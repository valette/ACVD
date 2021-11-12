/***************************************************************************
vtkCurvatureMeasure.h  -  description
-------------------
begin                : July 10 2005
copyright            : (C) 2005 by Sebastien Valette
email                : 
***************************************************************************/

/* ---------------------------------------------------------------------

* Copyright (c) CREATIS-LRMN (Centre de Recherche en Imagerie Medicale)
* Author : Sebastien Valette
*
*  This software is governed by the CeCILL-B license under French law and 
*  abiding by the rules of distribution of free software. You can  use, 
*  modify and/ or redistribute the software under the terms of the CeCILL-B 
*  license as circulated by CEA, CNRS and INRIA at the following URL 
*  http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html 
*  or in the file LICENSE.txt.
*
*  As a counterpart to the access to the source code and  rights to copy,
*  modify and redistribute granted by the license, users are provided only
*  with a limited warranty  and the software's author,  the holder of the
*  economic rights,  and the successive licensors  have only  limited
*  liability. 
*
*  The fact that you are presently reading this means that you have had
*  knowledge of the CeCILL-B license and that you accept its terms.
* ------------------------------------------------------------------------ */  

#include <vtkTimerLog.h>
#include <vtkTriangle.h>
#include <vtkCellData.h>
#include <vtkPriorityQueue.h>
#include <vtkMultiThreader.h>

#include "vtkCurvatureMeasure.h"
#include "vtkNeighbourhoodComputation.h"

#define DISPLAYINTERVAL 10000

/// this class is made to compute curvature measure of a vtkIdList with polynomial fitting
/// It was created to ease mumtithreading
class vtkSinglePolynomialMeasure:public vtkObject
{
public:

	void SetInputData (vtkSurface * Input)
	{
		this->Input = Input;
	};

	// Computes the curvature for a given group of triangles. Returns sqrt(c1^2+c2^2) where c1 and c2 are 
	// principal curvatures. Can also write the two principal directions and Curvatures if the second               parameter is provided. 
	// Principal directions and curvatures are stored this way : c1,x1,y1,z1,c2,x2,y2,z2;
	double ComputeFitting (vtkIdList * FList, double *CurvatureInfo = 0);

	vtkSurface *GetInput ()
	{
		return this->Input;
	};

	/// the public constructor
	static vtkSinglePolynomialMeasure *New ()
	{
		// First try to create the object from the vtkObjectFactory
		vtkObject *ret =
			vtkObjectFactory::
			CreateInstance ("vtkSinglePolynomialMeasure");
		if (ret)
		{
			return (vtkSinglePolynomialMeasure *) ret;
		}
		// If the factory was unable to create the object, then create it here.
		return (new vtkSinglePolynomialMeasure);
	}

	int NumberOfCellsWithSmallNeighbourhood;
	int NumberOfBadMatrices;

protected:

	vtkSinglePolynomialMeasure ()
	{
		this->NumberOfBadMatrices = 0;
		this->NumberOfCellsWithSmallNeighbourhood = 0;
		this->MatricesSize=0;
		this->AllocateMatrices(100);

		for (i = 0; i < 6; i++)
			Quadric[i] = new double;

		A = new double *[2];
		A[0] = new double[2];
		A[1] = new double[2];
		B = new double *[2];
		B[0] = new double[2];
		B[1] = new double[2];
		EigenValues = new double[2];
		EigenVectors = new double *[2];
		EigenVectors[0] = new double[2];
		EigenVectors[1] = new double[2];

	};

	~vtkSinglePolynomialMeasure ()
	{
		for (i = 0; i < 6; i++)
			delete Quadric[i];

		this->FreeMatrices();

		delete[]A[0];
		delete[]A[1];
		delete[]A;
		delete[]B[0];
		delete[]B[1];
		delete[]B;
		delete[]EigenVectors[0];
		delete[]EigenVectors[1];
		delete[]EigenVectors;
		delete[]EigenValues;
	};

private:

	/// The input vtkSurface;
	vtkSurface * Input;

	void AllocateMatrices(int Size)
	{
		this->FreeMatrices();

		SecondMember = new double * [Size];
		VandermondeMatrix = new double * [Size];
		for (i = 0; i < Size; i++)
		{
			VandermondeMatrix[i] = new double[6];
			SecondMember[i]= new double;
		}

		this->MatricesSize=Size;
	}

	void FreeMatrices()
	{
		if (this->MatricesSize==0)
			return;

		for (i = 0; i < this->MatricesSize; i++)
		{
			delete [] VandermondeMatrix[i];
			delete [] SecondMember[i];
		}

		delete [] VandermondeMatrix;
		delete [] SecondMember;
		this->MatricesSize=0;
	}

	int MatricesSize;

	double **VandermondeMatrix;
	double **SecondMember;
	double *Quadric[6];
	double **A, **B;
	double Barycenter[3];
	double *EigenValues, **EigenVectors;

	int i, j, k;
	double SArea, Area;
	int Face;
	vtkIdType v1, v2, v3;
	double Point1[3], Point2[3], Point3[3], Normal[3], Origin[3],
		Frame[3][3], x, y, z;

	// a version of vtkMath::SolveLeastsquares without memory leaks....
	int SolveLeastSquares (int numberOfSamples, double **xt, int xOrder,
			       double **yt, int yOrder, double **mt,
			       int checkHomogeneous);

};

double
vtkSinglePolynomialMeasure::ComputeFitting (vtkIdList * FList,
					    double *CurvatureInfo)
{
	SArea = 0;
	Origin[0] = 0;
	Origin[1] = 0;
	Origin[2] = 0;
	Frame[0][0] = 0;
	Frame[0][1] = 0;
	Frame[0][2] = 0;

	if (FList->GetNumberOfIds()>this->MatricesSize)
		this->AllocateMatrices(FList->GetNumberOfIds());

	// Compute mean normal, area and centroid of the region
	for (j = 0; j < FList->GetNumberOfIds (); j++)
	{
		Face = FList->GetId (j);
		this->GetInput ()->GetFaceVertices (Face, v1, v2, v3);
		this->GetInput ()->GetPointCoordinates (v1, Point1);
		this->GetInput ()->GetPointCoordinates (v2, Point2);
		this->GetInput ()->GetPointCoordinates (v3, Point3);
		vtkTriangle::ComputeNormal (Point1, Point2, Point3, Normal);
		this->Input->GetCellMassProperties (Face, Area, Barycenter);
		for (k = 0; k < 3; k++)
		{
			Origin[k] += Area * Barycenter[k];
			Frame[0][k] += Area * Normal[k];
		}
		SArea += Area;
	}

	// Normalization of the first basis vector

	for (j = 0; j < 3; j++)
	{
		Origin[j] = Origin[j] / SArea;
	}
	vtkMath::Normalize (Frame[0]);


	// Construction of the second and third basis vector for each frame

	Frame[2][0] = Frame[0][0];
	Frame[2][1] = Frame[0][1];
	Frame[2][2] = Frame[0][2];

	Frame[1][1] = Frame[2][0];
	Frame[1][2] = Frame[2][1];
	Frame[1][0] = Frame[2][2];

	vtkMath::Cross (Frame[1], Frame[2], Frame[0]);
	vtkMath::Normalize (Frame[0]);
	vtkMath::Cross (Frame[2], Frame[0], Frame[1]);

	if (FList->GetNumberOfIds () > 6)
	{

		// Construction of the Vandermonde matrix
		double h = 0;	// h is the order of magnitude of the coordinates x and y;

		for (j = 0; j < FList->GetNumberOfIds (); j++)
		{
			Face = FList->GetId (j);
			this->Input->GetCellMassProperties (Face, Area,
							    Barycenter);
			for (k = 0; k < 3; k++)
			{
				Barycenter[k] -= Origin[k];
			}

			x = vtkMath::Dot (Barycenter, Frame[0]);
			y = vtkMath::Dot (Barycenter, Frame[1]);
			z = vtkMath::Dot (Barycenter, Frame[2]);

			VandermondeMatrix[j][0] = 1.0;
			VandermondeMatrix[j][1] = x;
			VandermondeMatrix[j][2] = y;
			VandermondeMatrix[j][3] = x * x;
			VandermondeMatrix[j][4] = x * y;
			VandermondeMatrix[j][5] = y * y;

			SecondMember[j][0] = z;

			h += sqrt (x * x + y * y);
		}

		// Conditionning the Vandermonde Matrix

		h = h / ((double) FList->GetNumberOfIds ());


		for (j = 0; j < FList->GetNumberOfIds (); j++)
		{
			VandermondeMatrix[j][1] /= h;
			VandermondeMatrix[j][2] /= h;
			VandermondeMatrix[j][3] /= h * h;
			VandermondeMatrix[j][4] /= h * h;
			VandermondeMatrix[j][5] /= h * h;
		}
/*
#if ( (VTK_MAJOR_VERSION >= 5))
		if (vtkMath::
		    SolveLeastSquares (FList->GetNumberOfIds (),
				       VandermondeMatrix, 6, SecondMember, 1,
				       Quadric, 0) == 0)
#else
		if (vtkMath::
		    SolveLeastSquares (FList->GetNumberOfIds (),
				       VandermondeMatrix, 6, SecondMember, 1,
				       Quadric) == 0)
#endif
		*/

		if (this->SolveLeastSquares (FList->GetNumberOfIds (),
					     VandermondeMatrix, 6,
					     SecondMember, 1, Quadric,
					     0) == 0)
		{
			NumberOfBadMatrices++;
			if (CurvatureInfo)
			{
				// The curvature will be defined as 0 but we still need to give principal directions
				for (i = 0; i < 6; i++)
					CurvatureInfo[i] = 0;
			}

			return (0);
		}
		else
		{
			// invert the matrix conditionnement

			Quadric[1][0] /= h;
			Quadric[2][0] /= h;
			Quadric[3][0] /= h * h;
			Quadric[4][0] /= h * h;
			Quadric[5][0] /= h * h;

			double E, F, G, e, f, g;
			E = 1.0 + Quadric[1][0] * Quadric[1][0];
			F = Quadric[1][0] * Quadric[2][0];
			G = 1.0 + Quadric[2][0] * Quadric[2][0];
			e = 2.0 * Quadric[3][0] / sqrt (Quadric[1][0] *
							Quadric[1][0] + 1.0 +
							Quadric[2][0] *
							Quadric[2][0]);
			f = 2.0 * Quadric[4][0] / sqrt (Quadric[1][0] *
							Quadric[1][0] + 1.0 +
							Quadric[2][0] *
							Quadric[2][0]);
			g = 2.0 * Quadric[5][0] / sqrt (Quadric[1][0] *
							Quadric[1][0] + 1.0 +
							Quadric[2][0] *
							Quadric[2][0]);

			A[0][0] = E;
			A[0][1] = F;
			A[1][0] = F;
			A[1][1] = G;

			int test;
			test = vtkMath::InvertMatrix (A, B, 2);
			if (test == 0)
			{
				if (CurvatureInfo)
				{
					// The curvature will be defined as 0 but we still need to give principal directions
					for (i = 0; i < 6; i++)
						CurvatureInfo[i] = 0;

				}
				return (0);
			}

			// workaround for a good computation of curvature.... to be improved though
			if (1)
			{
				A[0][0] = -(e * B[0][0] + f * B[1][0]);
				A[1][0] = -(e * B[0][1] + f * B[1][1]);
				A[0][1] = -(f * B[0][0] + g * B[1][0]);
				A[1][1] = -(f * B[0][1] + g * B[1][1]);
			}
			else
			{
				A[0][0] = e ;
				A[1][0] = f;
				A[0][1] = f;
				A[1][1] = g ;
			}

			if (1)
			{
				test = 1;
				test = vtkMath::JacobiN (A, 2, EigenValues,
							 EigenVectors);
				if (test == 0)
				{
					if (CurvatureInfo)
					{
						// The curvature will be defined as 0 but we still need to give principal directions
						for (i = 0; i < 6; i++)
							CurvatureInfo[i] = 0;
					}
					return (0);
				}
			}
			else
			{

				// attempt to avoid the use of JacobiN (the matrix is actuelly not symetric)
				// not working correctly right now... if anyone has a solution, I buy it :)
				double A33[3][3], U[3][3], w[3], VT[3][3];

				for (i = 0; i < 3; i++)
				{
					for (j = 0; j < 3; j++)
					{
						if ((i < 2) && (j < 2))
							A33[i][j] = A[i][j];
						else
							A33[i][j] = 0;
					}
				}
				A33[2][2] = 1.0;

				vtkMath::SingularValueDecomposition3x3 (A33,
									U, w,
									VT);
				double MaxValue = 0;
				int IMax;

				for (i = 0; i < 3; i++)
				{
					if (fabs (MaxValue) < fabs (VT[i][2]))
					{
						IMax = i;
						MaxValue = VT[i][2];
					}
				}

				int Index = 0;

				for (i = 0; i < 3; i++)
				{
					if (i != IMax)
					{
						for (j = 0; j < 2; j++)
							EigenVectors[Index][j]
								= VT[i][j];
						EigenValues[Index] = w[i];
						Index++;
					}
				}
			}

			if (CurvatureInfo)
			{
				// We have to compute the princpial directions (from 2D to 3D)
				for (i = 0; i < 6; i++)
				{
					CurvatureInfo[i] = 0;
				}
				for (i = 0; i < 2; i++)
				{
					for (j = 0; j < 2; j++)
					{
						for (k = 0; k < 3; k++)
						{
							CurvatureInfo[k +
								      3 *
								      j] +=
								sqrt (fabs
								      (EigenValues
								       [j])) *
								EigenVectors
								[i][j] *
								Frame[i][k];
						}
					}
				}
				if (fabs (EigenValues[0]) <
				    fabs (EigenValues[1]))
				{
					double TempCurvature;
					for (i = 0; i < 3; i++)
					{
						TempCurvature =
							CurvatureInfo[i];
						CurvatureInfo[i] =
							CurvatureInfo[i + 3];
						CurvatureInfo[i + 3] =
							TempCurvature;
					}

				}

//                              vtkMath::Normalize (CurvatureInfo);
//                              vtkMath::Normalize (CurvatureInfo + 5);
			}
			return (sqrt
				(EigenValues[0] * EigenValues[0] +
				 EigenValues[1] * EigenValues[1]));
		}
	}
	else
	{
		NumberOfCellsWithSmallNeighbourhood++;
		if (CurvatureInfo)
		{
			// The curvature will be defined as 0 but we still need to give principal directions
			for (i = 0; i < 6; i++)
			{
				CurvatureInfo[i] = 0;
			}
		}

		return (0);
	}
}


int
vtkSinglePolynomialMeasure::SolveLeastSquares (int numberOfSamples,
					       double **xt, int xOrder,
					       double **yt, int yOrder,
					       double **mt,
					       int checkHomogeneous)
{
	// check dimensional consistency
	if ((numberOfSamples < xOrder) || (numberOfSamples < yOrder))
	{
		vtkGenericWarningMacro
			("Insufficient number of samples. Underdetermined.");
		return 0;
	}

	int i, j, k;

	// set up intermediate variables
	double **XXt = new double *[xOrder];	// size x by x
	double **XXtI = new double *[xOrder];	// size x by x
	double **XYt = new double *[xOrder];	// size x by y
	for (i = 0; i < xOrder; i++)
	{
		XXt[i] = new double[xOrder];
		XXtI[i] = new double[xOrder];

		for (j = 0; j < xOrder; j++)
		{
			XXt[i][j] = 0.0;
			XXtI[i][j] = 0.0;
		}

		XYt[i] = new double[yOrder];
		for (j = 0; j < yOrder; j++)
		{
			XYt[i][j] = 0.0;
		}
	}

	// first find the pseudoinverse matrix
	for (k = 0; k < numberOfSamples; k++)
	{
		for (i = 0; i < xOrder; i++)
		{
			// first calculate the XXt matrix, only do the upper half (symmetrical)
			for (j = i; j < xOrder; j++)
			{
				XXt[i][j] += xt[k][i] * xt[k][j];
			}

			// now calculate the XYt matrix
			for (j = 0; j < yOrder; j++)
			{
				XYt[i][j] += xt[k][i] * yt[k][j];
			}
		}
	}

	// now fill in the lower half of the XXt matrix
	for (i = 0; i < xOrder; i++)
	{
		for (j = 0; j < i; j++)
		{
			XXt[i][j] = XXt[j][i];
		}
	}

	// next get the inverse of XXt
	if (!(vtkMath::InvertMatrix (XXt, XXtI, xOrder)))
	{
		// clean up:
		// set up intermediate variables
		for (i = 0; i < xOrder; i++)
		{
			delete[]XXt[i];
			delete[]XXtI[i];

			delete[]XYt[i];
		}
		delete[]XXt;
		delete[]XXtI;
		delete[]XYt;
		return 0;
	}

	// next get m
	for (i = 0; i < xOrder; i++)
	{
		for (j = 0; j < yOrder; j++)
		{
			mt[i][j] = 0.0;
			for (k = 0; k < xOrder; k++)
			{
				mt[i][j] += XXtI[i][k] * XYt[k][j];
			}
		}
	}


	// clean up:
	// set up intermediate variables
	for (i = 0; i < xOrder; i++)
	{
		delete[]XXt[i];
		delete[]XXtI[i];

		delete[]XYt[i];
	}
	delete[]XXt;
	delete[]XXtI;
	delete[]XYt;
	return 1;
}

VTK_THREAD_RETURN_TYPE
vtkCurvatureMeasure::ThreadedCurvatureComputation (void *arg)
{
	vtkMultiThreader::ThreadInfo *Info = (vtkMultiThreader::ThreadInfo*) arg;

	vtkCurvatureMeasure *CurvatureMeasure = (vtkCurvatureMeasure *) Info->UserData;
	int MyId = Info->ThreadID;
	int NumberOfThreads=Info->NumberOfThreads;

	int Cell;
	int RingSize;
	double DistanceMax;
	int NeighbourhoodComputationMethod;
	vtkIdList *FList = vtkIdList::New ();

	int NumberOfCells;
	if (CurvatureMeasure->ElementsType == 0)
		NumberOfCells = CurvatureMeasure->Input->GetNumberOfCells ();
	else
		NumberOfCells = CurvatureMeasure->Input->GetNumberOfPoints ();


	vtkNeighbourhoodComputation *Neighbourhood =
		vtkNeighbourhoodComputation::New ();
	Neighbourhood->SetCellType (CurvatureMeasure->ElementsType);
	Neighbourhood->SetInputData (CurvatureMeasure->Input);

	vtkSinglePolynomialMeasure *Measure =
		vtkSinglePolynomialMeasure::New ();
	Measure->SetInputData (CurvatureMeasure->Input);

	RingSize = CurvatureMeasure->RingSize;
	DistanceMax = CurvatureMeasure->NeighbourhoodSize;
	NeighbourhoodComputationMethod = CurvatureMeasure->NeighbourhoodComputationMethod;
	double StartTime = CurvatureMeasure->StartTime;
	double ElapsedTime;
	double TotalTimeEstimated;


	for (Cell=MyId;Cell<NumberOfCells;Cell+=NumberOfThreads)
	{
		// Compute the neighbour cells of the cell i and stores the cells in the FList vtkIdList
		if (NeighbourhoodComputationMethod == 1)
		{
			Neighbourhood->ComputeDistanceRingCells (Cell,
								 DistanceMax,
								 FList);
		}
		else
		{
			// Compute the neighbourhood in the n-ring;
			Neighbourhood->ComputeNRingCells (Cell, RingSize,
							  FList);
		}

		if (CurvatureMeasure->CellsCurvatureInfo)
		{
			// We have to store the curvature infos (principal curvatures and directions)
			double CurvatureInfo[6];
			CurvatureMeasure->CellsCurvatureIndicator->SetValue (Cell,Measure->ComputeFitting(FList,CurvatureInfo));
			int i;
			for (i = 0; i < 6; i++)
				CurvatureMeasure->CellsCurvatureInfo->SetValue (6 *Cell+i,CurvatureInfo[i]);
		}
		else
			CurvatureMeasure->CellsCurvatureIndicator->SetValue (Cell,Measure->ComputeFitting(FList));


		if ((Cell % DISPLAYINTERVAL == 0) && (Cell != 0))
		{
			char CarriageReturn = 13;
#if ( (VTK_MAJOR_VERSION >= 5))
			ElapsedTime =CurvatureMeasure->Timer->GetUniversalTime () -StartTime;
#else
			ElapsedTime =CurvatureMeasure->Timer->GetCurrentTime () - StartTime;
#endif
			TotalTimeEstimated =ElapsedTime * NumberOfCells / Cell;
			cout << CarriageReturn << (int) (TotalTimeEstimated -ElapsedTime) <<
				" s remaining." << " Total time: " << (int)
				TotalTimeEstimated << " s. " << Cell /
				ElapsedTime << " Cells/second          " <<
				std::flush;
		}
	}
	CurvatureMeasure->StatisticsLock->lock ();
	CurvatureMeasure->NumberOfBadMatrices += Measure->NumberOfBadMatrices;
	CurvatureMeasure->NumberOfCellsWithSmallNeighbourhood +=
		Measure->NumberOfCellsWithSmallNeighbourhood;
	CurvatureMeasure->StatisticsLock->unlock ();
	FList->Delete ();
	Neighbourhood->Delete ();
	Measure->Delete ();
	return (VTK_THREAD_RETURN_VALUE);
}

vtkDataArrayCollection *
vtkCurvatureMeasure::GetCurvatureIndicator ()
{
	int i;
	int NumberOfElements;


	if (this->ElementsType == 0)
		NumberOfElements = this->Input->GetNumberOfCells ();
	else
		NumberOfElements = this->Input->GetNumberOfPoints ();

	if (!this->CellsCurvatureIndicator)
	{
		this->CellsCurvatureIndicator = vtkDoubleArray::New ();
		this->CellsCurvatureIndicator->SetNumberOfValues (NumberOfElements);
	}

	if (this->ComputeCurvatureInfoFlag)
	{
		if (!this->CellsCurvatureInfo)
		{
			this->CellsCurvatureInfo = vtkFloatArray::New ();
			this->CellsCurvatureInfo->SetNumberOfValues (NumberOfElements * 6);
		}
	}

	for (i = 0; i < NumberOfElements; i++)
		CellsCurvatureIndicator->SetValue (i, 0.0);


	vtkTimerLog *Timer = vtkTimerLog::New ();
	Timer->StartTimer ();

	cout << "Computing Curvature ........" << endl;

	if (this->ComputationMethod == 0)
		this->ComputeCurvatureIndicatorWithNormalAnalysis ();
	else
		this->ComputeCurvatureIndicatorWithPolynomialFitting ();


	Timer->StopTimer ();
	cout << "The curvature indicator calculation took :" << Timer->
		GetElapsedTime () << " seconds." << endl;
	Timer->Delete ();

	if (this->Display)
	{
		vtkDoubleArray *IndicatorColors = vtkDoubleArray::New ();
		IndicatorColors->SetNumberOfValues (NumberOfElements);
		double MaxIndicatorColor = 0;
		double MeanIndicator = 0;

		for (i = 0; i < NumberOfElements; i++)
		{
			IndicatorColors->SetValue (i,CellsCurvatureIndicator->GetValue (i));
			if (MaxIndicatorColor < IndicatorColors->GetValue (i))
				MaxIndicatorColor =	IndicatorColors->GetValue (i);
			MeanIndicator += IndicatorColors->GetValue (i);
		}
		MeanIndicator /= NumberOfElements;
		double ClampingFactor = 100.0;
		if (MaxIndicatorColor > ClampingFactor * MeanIndicator)
			MaxIndicatorColor = ClampingFactor * MeanIndicator;
		for (i = 0; i < NumberOfElements; i++)
		{
			if (IndicatorColors->GetValue (i) > MaxIndicatorColor)
				IndicatorColors->SetValue (i,MaxIndicatorColor);
			IndicatorColors->SetValue (i,IndicatorColors->GetValue(i)/MaxIndicatorColor);
		}

		RenderWindow *Window;
		vtkPolyData *Mesh2 = vtkPolyData::New ();
		Window = RenderWindow::New ();

		Mesh2->ShallowCopy (this->Input);
		if (this->ElementsType == 0)
		{
			Mesh2->GetCellData ()->SetScalars (IndicatorColors);
			Mesh2->GetPointData ()->SetScalars (0);
		}
		else
		{
			Mesh2->GetPointData ()->SetScalars (IndicatorColors);
			Mesh2->GetCellData ()->SetScalars (0);
		}

		Window->SetInputData (Mesh2);

		vtkLookupTable *bwLut = vtkLookupTable::New ();
		bwLut->SetTableRange (0, 1);

		bwLut->SetSaturationRange (0, 0);	// Black     and     white
		bwLut->SetHueRange (0, 0);
		bwLut->SetValueRange (0, 1);
		bwLut->SetScaleToLog10 ();
		bwLut->Build ();
		Window->SetLookupTable (bwLut);

		Window->Render ();
		Window->SetWindowName ("Curvature Measure");
		if (this->AnchorRenderWindow)
			Window->AttachToRenderWindow (AnchorRenderWindow);
		else
			this->AnchorRenderWindow = Window;
		Window->Interact ();

		bwLut->Delete ();
	}

	this->CurvatureCollection = vtkDataArrayCollection::New ();
	this->CurvatureCollection->AddItem (this->CellsCurvatureIndicator);
	if (this->CellsCurvatureInfo)
		this->CurvatureCollection->AddItem (this->CellsCurvatureInfo);

	return (this->CurvatureCollection);
}

vtkPolyData *
vtkCurvatureMeasure::GetPrincipalDirectionsPolyData ()
{

	int i;
	vtkPolyData *Lines = vtkPolyData::New ();
	vtkIntArray *Colors = vtkIntArray::New ();


	int NumberOfElements;

	if (this->ElementsType == 0)
	{
		NumberOfElements = this->Input->GetNumberOfCells ();
	}
	else
	{
		NumberOfElements = this->Input->GetNumberOfPoints ();

	}
	Colors->SetNumberOfValues (2 * NumberOfElements);

	Lines->Allocate (NumberOfElements * 2);
	vtkPoints *Points = vtkPoints::New ();
	Points->Allocate (NumberOfElements * 4);
	Lines->SetPoints (Points);
	Points->Delete ();

	vtkDoubleArray *EdgesLength = this->Input->GetEdgeLengths ();

	double AverageLength = 0;
	for (i = 0; i < this->Input->GetNumberOfEdges (); i++)
	{
		AverageLength += EdgesLength->GetValue (i);
	}
	AverageLength /= (double) this->Input->GetNumberOfEdges ();


	vtkIdList *VList = vtkIdList::New ();
	VList->SetNumberOfIds (2);
	for (i = 0; i < NumberOfElements; i++)
	{
		double Center[3];
		this->Input->GetPoint (i, Center);

		double P1[3], P2[3], P3[3], P4[3];

		float *CurvatureInfo;
		CurvatureInfo = this->CellsCurvatureInfo->GetPointer (i * 6);
		int j;

		double C1, C2;
		C1 = 0.1;	//(fabs (vtkMath::Norm(CurvatureInfo)));
		C2 = 0.1;	//(fabs (vtkMath::Norm(CurvatureInfo+3)));


		double Factor = AverageLength;

		for (j = 0; j < 3; j++)
		{
			P1[j] = Center[j] - CurvatureInfo[j] * C1 * Factor;
			P2[j] = Center[j] + CurvatureInfo[j] * C1 * Factor;
			P3[j] = Center[j] - CurvatureInfo[j +
							  3] * C2 * Factor;
			P4[j] = Center[j] + CurvatureInfo[j +
							  3] * C2 * Factor;
		}
		int v1, v2, v3, v4;
		v1 = Lines->GetPoints ()->InsertNextPoint (P1);
		v2 = Lines->GetPoints ()->InsertNextPoint (P2);
		v3 = Lines->GetPoints ()->InsertNextPoint (P3);
		v4 = Lines->GetPoints ()->InsertNextPoint (P4);

		VList->SetId (0, v1);
		VList->SetId (1, v2);
		Lines->InsertNextCell (VTK_LINE, VList);
		Colors->SetValue (2 * i, 0);
		Colors->SetValue (2 * i + 1, 1);

		VList->SetId (0, v3);
		VList->SetId (1, v4);
		Lines->InsertNextCell (VTK_LINE, VList);
	}
	Lines->GetCellData ()->SetScalars (Colors);
	VList->Delete();
	return (Lines);
}


void
vtkCurvatureMeasure::ComputeCurvatureIndicatorWithPolynomialFitting ()
{
	double Point1[3];
	double Bounds[6];
	this->GetInput ()->GetPoints ()->ComputeBounds ();
	this->GetInput ()->GetPoints ()->GetBounds (Bounds);

	Point1[0] = Bounds[1] - Bounds[0];
	Point1[1] = Bounds[3] - Bounds[2];
	Point1[2] = Bounds[5] - Bounds[4];

	double DistanceMax;
	DistanceMax = vtkMath::Norm (Point1) * this->NeighbourhoodSize;

	int NumberOfElements;
	if (this->ElementsType == 0)
		NumberOfElements = this->Input->GetNumberOfCells ();
	else
		NumberOfElements = this->Input->GetNumberOfPoints ();

#if ( (VTK_MAJOR_VERSION >= 5))
	this->StartTime = this->Timer->GetUniversalTime ();
#else
	this->StartTime = this->Timer->GetCurrentTime ();
#endif
	vtkMultiThreader *Threader=vtkMultiThreader::New();
		
	Threader->SetSingleMethod (ThreadedCurvatureComputation, (void *) this);
	if (this->NumberOfThreads)
		Threader->SetNumberOfThreads (this->NumberOfThreads);
	Threader->SingleMethodExecute ();
	cout << (char) 13;
	
	Threader->Delete();
			
	if (NumberOfBadMatrices > 0)
	{
		cout << endl << NumberOfBadMatrices <<
			" bad matrices for the computation of the curvature"
			<< endl;
	}

	if (NumberOfCellsWithSmallNeighbourhood > 0)
		cout << endl<< NumberOfCellsWithSmallNeighbourhood <<
			" cells with too small neighbourhood." << endl;

}

void
vtkCurvatureMeasure::ComputeCurvatureIndicatorWithNormalAnalysis ()
{
	vtkPriorityQueue *EigenValuesQueue = vtkPriorityQueue::New ();
	EigenValuesQueue->Allocate (3);


	vtkIdType i = 0, j, k, l;
	vtkIdType v1, v2, v3, Edge1, f1, f2;
	vtkIdList *VList = vtkIdList::New ();
	vtkIdList *VList2 = vtkIdList::New ();
	vtkIdList *FList = vtkIdList::New ();
	vtkIdList *EList = vtkIdList::New ();
	vtkIdList *FList2 = vtkIdList::New ();
	vtkIdList *EList2 = vtkIdList::New ();
	vtkIdList *TempList;
	double Area, SArea;
	double Tensor[3][3], EigenVectors[3][3], EigenValues[3];
	double Temp[3];


	vtkIntArray *VisitedCells = vtkIntArray::New ();
	VisitedCells->SetNumberOfValues (this->GetInput ()->
					 GetNumberOfCells ());
	for (i = 0; i < this->GetInput ()->GetNumberOfCells (); i++)
		VisitedCells->SetValue (i, -1);

	vtkIntArray *VisitedVertices = vtkIntArray::New ();
	VisitedVertices->SetNumberOfValues (this->Input->
					    GetNumberOfPoints ());
	for (i = 0; i < this->GetInput ()->GetNumberOfPoints (); i++)
		VisitedVertices->SetValue (i, -1);

	for (i = 0; i < this->Input->GetNumberOfCells (); i++)
	{
		this->Input->GetFaceVertices (i, v1, v2, v3);
		VList->Reset ();
		VList->InsertNextId (v1);
		VList->InsertNextId (v2);
		VList->InsertNextId (v3);

		FList->Reset ();
		FList->InsertNextId (i);
		VisitedCells->SetValue (i, i);
		for (j = 0; j < RingSize; j++)
		{
			VList2->Reset ();
			for (k = 0; k < VList->GetNumberOfIds (); k++)
			{
				v1 = VList->GetId (k);
				if (VisitedVertices->GetValue (v1) != i)
				{
					VisitedVertices->SetValue (v1, i);
					this->Input->
						GetVertexNeighbourEdges (v1,
									 EList);
					for (l = 0;
					     l < EList->GetNumberOfIds ();
					     l++)
					{
						Edge1 = EList->GetId (l);
						this->Input->
							GetEdgeVertices
							(Edge1, v2, v3);
						if (v1 == v2)
							VList2->InsertNextId
								(v3);
						else
							VList2->InsertNextId
								(v2);

						this->Input->
							GetEdgeFaces (Edge1,
								      f1, f2);
						if (VisitedCells->
						    GetValue (f1) != i)
						{
							FList->InsertNextId
								(f1);
							VisitedCells->
								SetValue (f1,
									  i);
						}
						if (f2 >= 0)
						{
							if (VisitedCells->
							    GetValue (f2) !=
							    i)
							{
								FList->InsertNextId (f2);
								VisitedCells->
									SetValue
									(f2,
									 i);
							}
						}
					}
				}

			}
			TempList = VList;
			VList = VList2;
			VList2 = TempList;
		}


		// Compute the tensor

		for (j = 0; j < 3; j++)
		{
			for (k = 0; k < 3; k++)
			{
				Tensor[j][k] = 0;
			}
		}
		// Measure Area of the part
		SArea = 0;
		vtkIdType v1, v2, v3;
		double P1[3], P2[3], P3[3], N[3];

		for (j = 0; j < FList->GetNumberOfIds (); j++)
		{
			this->Input->GetCellMassProperties (FList->GetId (j),
							    Area, Temp);

			this->GetInput ()->GetFaceVertices (FList->GetId (j),
							    v1, v2, v3);


			this->GetInput ()->GetPointCoordinates (v1, P1);
			this->GetInput ()->GetPointCoordinates (v2, P2);
			this->GetInput ()->GetPointCoordinates (v3, P3);

			vtkTriangle::ComputeNormal (P1, P2, P3, N);


			SArea += Area;
			for (k = 0; k < 3; k++)
			{
				for (l = 0; l < 3; l++)
				{
					Tensor[k][l] += Area * N[k] * N[l];
				}
			}
		}

		// divide the tensor by the area;
		if (SArea)
		{
			for (j = 0; j < 3; j++)
			{
				for (k = 0; k < 3; k++)
				{
					Tensor[j][k] = Tensor[j][k] / SArea;
				}
			}
		}

		vtkMath::Diagonalize3x3 (Tensor, EigenValues, EigenVectors);

		// Sort the eigen values and eigen vectors (descendant ranking)
		EigenValuesQueue->Reset ();
		for (j = 0; j < 3; j++)
		{
			EigenValuesQueue->Insert (-fabs (EigenValues[j]), j);
		}
		v1 = EigenValuesQueue->Pop ();
		v2 = EigenValuesQueue->Pop ();
		v3 = EigenValuesQueue->Pop ();
		CellsCurvatureIndicator->SetValue (i,sqrt (EigenValues[v2]*EigenValues[v2]+EigenValues[v3]*EigenValues[v3]));

	}

	VList->Delete ();
	VList2->Delete ();
	FList->Delete ();
	EList->Delete ();
	FList2->Delete ();
	EList2->Delete ();
	VisitedCells->Delete ();
	VisitedVertices->Delete ();

}

vtkCurvatureMeasure *
vtkCurvatureMeasure::New ()
{
	// First try to create the object from the vtkObjectFactory
	vtkObject *ret =
		vtkObjectFactory::CreateInstance ("vtkCurvatureMeasure");
	if (ret)
	{
		return (vtkCurvatureMeasure *) ret;
	}
	// If the factory was unable to create the object, then create it here.
	return (new vtkCurvatureMeasure);
}

vtkCurvatureMeasure::vtkCurvatureMeasure ()
{
	this->CellsCurvatureIndicator = 0;
	this->AnchorRenderWindow = 0;
	this->ComputationMethod = 1;
	this->NeighbourhoodComputationMethod = 0;
	this->RingSize = 3;
	this->NeighbourhoodSize = 0.01;
	this->NumberOfThreads = 0;
	this->Display = 0;

	this->ElementsType = 0;

	this->DisplayCurvatureInfoFlag = 0;
	this->ComputeCurvatureInfoFlag = 1;
	this->CellsCurvatureInfo = 0;
	
	this->StatisticsLock = new std::mutex();
	this->Timer = vtkTimerLog::New ();
	this->NumberOfBadMatrices = 0;
	this->NumberOfCellsWithSmallNeighbourhood = 0;	
}

vtkCurvatureMeasure::~vtkCurvatureMeasure ()
{
	delete this->StatisticsLock;
	this->Timer->Delete ();
	
	if (this->CurvatureCollection)
		this->CurvatureCollection->Delete();
	if (this->CellsCurvatureIndicator)
		this->CellsCurvatureIndicator->Delete();
	if (this->CellsCurvatureInfo)
		this->CellsCurvatureInfo->Delete();
}
