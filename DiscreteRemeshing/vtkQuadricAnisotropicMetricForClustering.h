/*=========================================================================

  Program:   vtkAnisotropicMetricForClustering.h
  Language:  C++
  Date:      2006/03
  Auteur:    Sebastien VALETTE

=========================================================================*/

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

#ifndef _VTKQUADRICANISOTROPICMETRICFORCLUSTERING_H_
#define _VTKQUADRICANISOTROPICMETRICFORCLUSTERING_H_

#include <vtkTriangle.h>
#include <vtkMath.h>

#include <vtkTriangle.h>
#include <vtkFloatArray.h>
#include "vtkSurface.h"
#include "vtkQuadricTools.h"
#define	AnisotropicSVTHRESHOLD	0.005

class vtkQuadricAnisotropicMetricForClustering
{
public:

	void MultiplyItemWeight(vtkIdType ItemId, double Factor)
	{
		for (int i=0;i<3;i++)
			this->Items[ItemId].Value[i]*=Factor;
		this->Items[ItemId].Weight*=Factor;
	}

	void SetQuadricsOptimizationLevel(int L)
	{
		this->QuadricsOptimizationLevel=L;
	}

	// to be implemented
	void SetConstrainedClustering(int C){}

	int IsCurvatureIndicatorNeeded ()
	{
		return 1;
	}
	int IsPrincipalDirectionsNeeded ()
	{
		return (1);
	}

	void SetCurvatureInfo (vtkDataArrayCollection * Info)
	{
		if (this->CustomWeights)
			this->CustomWeights->Delete();
		this->CustomWeights=(vtkDoubleArray*) Info->GetItem(0);
		this->CustomWeights->Register(this->Object);
		
		if (this->PrincipalDirections)
			this->PrincipalDirections->Delete();
		this->PrincipalDirections=(vtkFloatArray*) Info->GetItem(1);
		this->PrincipalDirections->Register(this->Object);
	}
	void SetGradation (double Gradation)
	{
		this->Gradation = Gradation;
	}

	double GetGradation ()
	{
		return this->Gradation;
	}

	struct Item
	{
		float Value[3];
		float Weight;

		// the 3x3 matrix (in compact form because it is symetric)
		// of the Riemmanian tensor used to compute distances.
		// Values are stored this way: [T11 T12 T13 T22 T23 T33]
		// this tensor is symetric;
		float *Tensor;

		float TensorXCentroid[3];

		// the 4x4 QEM matrix in compact form
		// note : there are only 9 coefficients as the last one is useless
		// in this context
		double Quadric[9];

	};

	void ResetItem (Item * I)
	{
		I->Value[0] = 0;
		I->Value[1] = 0;
		I->Value[2] = 0;
		I->Weight = 0;
		int i;
		for (i = 0; i < 6; i++)
			I->Tensor[i] = 0;

		for (i = 0; i < 3; i++)
			I->TensorXCentroid[i] = 0;

		for (i = 0; i < 9; i++)
			I->Quadric[i] = 0;
	};

	double GetItemWeight (vtkIdType ItemId)
	{
		return this->Items[ItemId].Weight;
	};

	struct Cluster
	{
		double SValue[3];
		double SWeight;
		double EnergyValue;

		// the 3x3 Accumulation matrix (in compact form because it is symetric)
		// of the Riemmanian tensors used to compute distances.
		// Values are stored this way: [ST11 ST12 ST13 ST22 ST23 ST33]
		// this tensor is symetric;
		double STensor[6];

		// Accumulation vector for the sum of Tensors Multiplied by the
		// Centroids
		double STensorXCentroid[3];

		double Centroid[3];

		// The accumulated quadrics
		double SQuadric[9];

	};

	int GetClusterRankDeficiency (Cluster * C)
	{
		return (0);
	};

	double GetClusterEnergy (Cluster * C)
	{
		return C->EnergyValue;
	}
	void ComputeClusterEnergy (Cluster * C)
	{
		double x, y, z;
		x = C->Centroid[0];
		y = C->Centroid[1];
		z = C->Centroid[2];

		double *T = C->STensor;

		C->EnergyValue = T[0] * x * x + T[3] * y * y + T[5] * z * z
			+ 2.0 * T[1] * x * y + 2.0 * T[2] * x * z +
			2.0 * T[4] * y * z;

		T = C->STensorXCentroid;

		C->EnergyValue -= 2.0 * (x * T[0] + y * T[1] + z * T[2]);
	}
	double ComputeDistanceBetweenItemAndCluster (Item * I, Cluster * C)
	{
		int i;
		double Cluster[3];

		for (i = 0; i < 3; i++)
		{
			Cluster[i] = C->SValue[i] / C->SWeight;
		}
		return (this->
			ComputeDistanceBetweenItemAndPoint (I, Cluster));
	}

	double ComputeDistanceBetweenItemAndPoint (Item * I, double *P1)
	{
		cout<<"Not Implemented!"<<endl;

		return (0);
	}


	void DeepCopy (Cluster * Source, Cluster * Destination)
	{
		int i;
		Destination->SWeight = Source->SWeight;
		for (i = 0; i < 3; i++)
		{
			Destination->SValue[i] = Source->SValue[i];
			Destination->STensorXCentroid[i] =
				Source->STensorXCentroid[i];
			Destination->Centroid[i]=Source->Centroid[i];
		}

		for (i = 0; i < 6; i++)
			Destination->STensor[i] = Source->STensor[i];

		for (i = 0; i < 9; i++)
			Destination->SQuadric[i] = Source->SQuadric[i];
		Destination->EnergyValue = Source->EnergyValue;
	}

	void Add(Cluster *Source, vtkIdType ItemId, Cluster *Destination)
	{
		this->DeepCopy(Source,Destination);
		this->AddItemToCluster(ItemId, Destination);
	}

	void AddItemToCluster (vtkIdType ItemId, Cluster * C)
	{
		Item *I=this->Items+ItemId;
		int i;
		for (i = 0; i < 3; i++)
		{
			C->SValue[i] +=(double) I->Value[i];
			C->STensorXCentroid[i] +=(double) I->TensorXCentroid[i];
		}
		for (i = 0; i < 6; i++)
			C->STensor[i] +=(double) I->Tensor[i];

		for (i = 0; i < 9; i++)
			C->SQuadric[i] +=(double) I->Quadric[i];

		C->SWeight +=(double) I->Weight;
	}

	void Sub(Cluster *Source, vtkIdType ItemId, Cluster *Destination)
	{
		this->DeepCopy(Source,Destination);
		this->SubstractItemFromCluster(ItemId, Destination);		
	}

	void SubstractItemFromCluster (vtkIdType ItemId, Cluster * C)
	{
		Item *I=this->Items+ItemId;
		int i;
		for (i = 0; i < 3; i++)
		{
			C->SValue[i] -=(double) I->Value[i];
			C->STensorXCentroid[i] -=(double) I->TensorXCentroid[i];
		}
		for (i = 0; i < 6; i++)
			C->STensor[i] -=(double) I->Tensor[i];

		for (i = 0; i < 9; i++)
			C->SQuadric[i] -=(double) I->Quadric[i];

		C->SWeight -=(double) I->Weight;

	}

	void GetClusterCentroid (Cluster * C, double *P)
	{
		P[0] = C->Centroid[0];
		P[1] = C->Centroid[1];
		P[2] = C->Centroid[2];
	}

	void ComputeClusterCentroid (Cluster * C)
	{		
		C->Centroid[0]=C->SValue[0]/C->SWeight;
		C->Centroid[1]=C->SValue[1]/C->SWeight;
		C->Centroid[2]=C->SValue[2]/C->SWeight;
		vtkQuadricTools::ComputeRepresentativePoint(C->SQuadric,C->Centroid,this->QuadricsOptimizationLevel);
	}

	void SetClusterCentroid (Cluster * C, double *P)
	{
		C->SValue[0] = P[0] * C->SWeight;
		C->SValue[1] = P[1] * C->SWeight;
		C->SValue[2] = P[2] * C->SWeight;
	}
	void ResetCluster (Cluster * C)
	{
		int i;
		for (i = 0; i < 3; i++)
		{
			C->SValue[i] = 0;
			C->STensorXCentroid[i] = 0;
			C->Centroid[i]=0;
		}
		C->SWeight = 0;
		C->EnergyValue = 0;

		for (i = 0; i < 6; i++)
			C->STensor[i] = 0;

		for (i = 0; i < 9; i++)
			C->SQuadric[i] = 0;
	}
	vtkQuadricAnisotropicMetricForClustering ()
	{
		this->CustomWeights = 0;
		this->PrincipalDirections = 0;
		this->Gradation = 0;
		this->Object=vtkObject::New();
		this->QuadricsOptimizationLevel=3;
	};
	~vtkQuadricAnisotropicMetricForClustering ()
	{
		if (this->CustomWeights)
			this->CustomWeights->Delete();
		if (this->PrincipalDirections)
			this->PrincipalDirections->Delete();
		this->Object->Delete();	
	};

	void BuildMetric ( vtkSurface * Mesh, int ClusteringType ) {

		vtkIdType i;
		int NumberOfElements;

		//Build the items

		if (ClusteringType == 0)
		{
			// Items are triangles
			vtkIdType v1, v2, v3;
			double P1[3], P2[3], P3[3];
			NumberOfElements = Mesh->GetNumberOfCells ();
			Items = new Item[Mesh->GetNumberOfCells ()];
			for (i = 0; i < Mesh->GetNumberOfCells (); i++)
			{
				this->ComputeTriangleQuadric(Items+i,Mesh,i);
				
				Mesh->GetFaceVertices (i, v1, v2, v3);
				Mesh->GetPoint (v1, P1);
				Mesh->GetPoint (v2, P2);
				Mesh->GetPoint (v3, P3);
				if (Gradation == 0)
					Items[i].Weight =
						vtkTriangle::
						TriangleArea (P1, P2, P3);
				else
					Items[i].Weight =
						vtkTriangle::
						TriangleArea (P1, P2,
							      P3) *
						pow (CustomWeights->
						     GetValue (i), Gradation);
				Items[i].Value[0] =
					(P1[0] + P2[0] + P3[0]) / 3.0;
				Items[i].Value[1] =
					(P1[1] + P2[1] + P3[1]) / 3.0;
				Items[i].Value[2] =
					(P1[2] + P2[2] + P3[2]) / 3.0;
			}
		}
		else
		{
			// Items are vertices
			NumberOfElements = Mesh->GetNumberOfPoints ();
			
			double Point[3];
			

			Items = new Item[Mesh->GetNumberOfPoints ()];
			for (i = 0; i < Mesh->GetNumberOfPoints (); i++)
			{
				this->ComputeVertexQuadric(Items+i,Mesh,i);
				if (Gradation == 0)
					Items[i].Weight =
						Mesh->GetVertexArea (i);
				else
					Items[i].Weight =
						Mesh->
						GetVertexArea (i) *
						pow (CustomWeights->
						     GetValue (i), Gradation);

				Mesh->GetPoint (i, Point);
				Items[i].Value[0]=Point[0];
				Items[i].Value[1]=Point[1];
				Items[i].Value[2]=Point[2];
			}
		}

		this->ClampWeights (Items, NumberOfElements, 100000);

		for (i = 0; i < NumberOfElements; i++)
		{

			// Compute the tensor in compact form
			double Area;

			if (ClusteringType == 0)
				Area = Mesh->GetFaceArea (i);
			else
				Area = Mesh->GetVertexArea (i);

			double PrincipalDirections[6];
			int j;
			for (j = 0; j < 6; j++)
				PrincipalDirections[j] =
					this->PrincipalDirections->
					GetValue (i * 6 + j);

			Items[i].Tensor =
				this->PrincipalDirections->GetPointer (i * 6);

			Items[i].Tensor[0] =
				Area * PrincipalDirections[0] *
				PrincipalDirections[0] +
				Area * PrincipalDirections[3] *
				PrincipalDirections[3];

			Items[i].Tensor[1] =
				Area * PrincipalDirections[0] *
				PrincipalDirections[1] +
				Area * PrincipalDirections[3] *
				PrincipalDirections[4];

			Items[i].Tensor[2] =
				Area * PrincipalDirections[0] *
				PrincipalDirections[2] +
				Area * PrincipalDirections[3] *
				PrincipalDirections[5];

			Items[i].Tensor[3] =
				Area * PrincipalDirections[1] *
				PrincipalDirections[1] +
				Area * PrincipalDirections[5] *
				PrincipalDirections[5];

			Items[i].Tensor[4] =
				Area * PrincipalDirections[1] *
				PrincipalDirections[2] +
				Area * PrincipalDirections[4] *
				PrincipalDirections[5];

			Items[i].Tensor[5] =
				Area * PrincipalDirections[2] *
				PrincipalDirections[2] +
				Area * PrincipalDirections[5] *
				PrincipalDirections[5];


			// Compute the TensorXCentroid vector

			Items[i].TensorXCentroid[0] =
				Items[i].Tensor[0] * Items[i].Value[0]
				+
				Items[i].Tensor[1] *
				Items[i].Value[1] +
				Items[i].Tensor[2] * Items[i].Value[2];

			Items[i].TensorXCentroid[1] =
				Items[i].Tensor[1] * Items[i].Value[0]
				+
				Items[i].Tensor[3] *
				Items[i].Value[1] +
				Items[i].Tensor[4] * Items[i].Value[2];

			Items[i].TensorXCentroid[2] =
				Items[i].Tensor[2] * Items[i].Value[0]
				+
				Items[i].Tensor[4] *
				Items[i].Value[1] +
				Items[i].Tensor[5] * Items[i].Value[2];

			Items[i].Value[0] *= Items[i].Weight;
			Items[i].Value[1] *= Items[i].Weight;
			Items[i].Value[2] *= Items[i].Weight;
		}
	}
	// this method clamps the weights between AverageValue/Ratio and AverageValue*Ratio
	static void ClampWeights (Item * Items, int NumberOfValues,
				double Ratio)
	{
		double Average = 0;
		int i;
		for (i = 0; i < NumberOfValues; i++)
		{
			Average += Items[i].Weight;
		}
		Average = Average / (double) NumberOfValues;
		double Min = Average / Ratio;
		double Max = Average * Ratio;
		for (i = 0; i < NumberOfValues; i++)
		{
			if (Items[i].Weight > Max)
				Items[i].Weight = Max;
			if (Items[i].Weight < Min)
				Items[i].Weight = Min;
		}
	}

private:
	int QuadricsOptimizationLevel;
	vtkDoubleArray * CustomWeights;
	vtkFloatArray *PrincipalDirections;
	double Gradation;

	// Dummy object used for registering the curvature indicators
	vtkObject *Object;
	
	// the array storing the items
	Item *Items;	
	

	void ComputeVertexQuadric (Item * I, vtkSurface * Mesh, int Vertex)
	{
		
		int i;

		for (i=0;i<9;i++)
			I->Quadric[i]=0;

		vtkIdList *FList = vtkIdList::New ();
		Mesh->GetVertexNeighbourFaces (Vertex, FList);
		for (i = 0; i < FList->GetNumberOfIds (); i++)
		{
			vtkQuadricTools::AddTriangleQuadric(I->Quadric,Mesh,FList->GetId(i),false);                
		}
		FList->Delete ();
	}
	

	void ComputeTriangleQuadric (Item * I, vtkSurface * Mesh, int Face)
	{
		int i;
		for (i=0;i<9;i++)
			I->Quadric[i]=0;
		vtkQuadricTools::AddTriangleQuadric(I->Quadric,Mesh,Face,false);
	}
};


#endif
