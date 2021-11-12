/***************************************************************************
                          vtkCurvatureMeasure.h  -  description
                             -------------------
    begin                : July 10 2005
    copyright            : (C) 2003 by Sebastien Valette
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

#ifndef _VTKCURVATUREMEASURE_H_
#define _VTKCURVATUREMEASURE_H_

#include <vtkObjectFactory.h>
#include <vtkDataArrayCollection.h>
#include <vtkFloatArray.h>
#include <mutex>
#include "vtkSurface.h"
#include "RenderWindow.h"
class VTK_EXPORT vtkCurvatureMeasure : public vtkObject
{
public:

	/// the public constructor
	static vtkCurvatureMeasure* New();

	/// To define the input Mesh
	void SetInputData(vtkSurface *Input) {this->Input=Input;};

	/// Returns the Input Mesh
	vtkSurface *GetInput() {return this->Input;};

	/// Computes the curvature measure. 
	/// The results are embedded in a vtkDataArrayCollection.
	/// The first vtkDataArray is a double array containing the curvature indicator
	/// The (possibly) second array is a double array containing principal directions and curvatures
	vtkDataArrayCollection *GetCurvatureIndicator();
	
	// Returns a vtkPolyData containing the lines representing the principal directions
	vtkPolyData *GetPrincipalDirectionsPolyData();
	
	
	// Defines whether the principal directions and curvatures will be 
	// stored in the output vtkDataArrayCollection
	void SetComputePrincipalDirections(int S) {this->ComputeCurvatureInfoFlag=S;};

	// Set/ get the type of elements to compute curvature on (0 : triangles; 1 : vertices)
	vtkSetMacro(ElementsType,int);
	vtkGetMacro(ElementsType,int);

	vtkSetMacro(ComputationMethod,int);
	vtkGetMacro(ComputationMethod,int);

	vtkSetMacro(NeighbourhoodComputationMethod,int);
	vtkGetMacro(NeighbourhoodComputationMethod,int);
	
	// defines the size of the ring used in neighbourhood computation
	vtkSetMacro(RingSize,int);
	vtkGetMacro(RingSize,int);

	vtkSetMacro(NeighbourhoodSize,double);
	vtkGetMacro(NeighbourhoodSize,double);
	
	void SetDisplay(int D) {this->Display=D;};
	void SetNumberOfThreads(int N) {this->NumberOfThreads=N;};

protected:

	/// Timer for performances measures
	vtkTimerLog *Timer;
	double StartTime;

	// This Mutex is used for the gathering of statistics for the curvature measure
	// Statistics are the number of matrices with bad conditionment and the number of cells with too small neighbourhood
	std::mutex *StatisticsLock;
	int NumberOfBadMatrices;
	int NumberOfCellsWithSmallNeighbourhood;

	vtkCurvatureMeasure();
	~vtkCurvatureMeasure();

private:

	// the type of elements to compute curvature on.
	// 0: faces 1: vertices
	int ElementsType;

	// Value storing the curvature measure computation method.
	// 0: Normal analysis. 1:Polynomial fitting
	int ComputationMethod;
	
	// Flag defining wether the curvatures info will be stored instead of the simple curvature indicator
	int ComputeCurvatureInfoFlag;
	
	// Flag defining whether the principal directions will be displayed or not
	int DisplayCurvatureInfoFlag;

	// Value storing the neighbourhood coputation for each cell.
	// 0: Compute the n-ring.  1: Compute the neighbourhood (according to the Bounding Box Diagonal).
	int NeighbourhoodComputationMethod;

	// Value storing the Ringsize when NeighbourhoodComputationMethod=0.
	int RingSize;

	// Value storing the Neighbourhood size when NeighbourhoodComputationMethod=1.
	double NeighbourhoodSize;

	// The principal directions
	vtkFloatArray *CellsCurvatureInfo;
	
	// The curvature indicator :  sqrt (c1*c1+c2*c2)
	vtkDoubleArray *CellsCurvatureIndicator;

	// The Collection containing both Curvature indicator and principal directions	
	vtkDataArrayCollection *CurvatureCollection;
	
	RenderWindow *AnchorRenderWindow;

	/// this method computes a curvature indicator having this formula : sqrt(c1*c1+c2*c2)
	/// where c1 and c2 are the local curvature measures, computed by polynomial fitting, as explained in:
	/// Estimating Differential Quantities using Polynomial fitting of Osculating Jets", Cazals and Pouget, 2005
	void ComputeCurvatureIndicatorWithPolynomialFitting();

	/// this method computes a curvature indicator based on the proposed solution in:
	/// "Adaptive polygonal mesh simplification with discrete centroidal Voronoi Tesselations", Valette, Kompastsiaris, Chassery, 2005
	void ComputeCurvatureIndicatorWithNormalAnalysis();

	/// The input vtkSurface;
	vtkSurface *Input;

	// the number of threads used for the computation (used only when computing polynomial fitting)
	// Default value is set to the number of processors
	int NumberOfThreads;

	// parameter defining whether the curvature indicator will be displayed or not
	int Display;
	
	// the threaded method to compute the curvature
	static VTK_THREAD_RETURN_TYPE ThreadedCurvatureComputation (void *arg);

};

#endif
