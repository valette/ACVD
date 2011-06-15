/*=========================================================================

  Program:   vtkManifoldSurfaceSimplification
  Module:    vtkSurface
  Language:  C++
  Date:      2008/03
  Auteur:    Sebastien VALETTE
  
=========================================================================*/

#ifndef __vtkManifoldSimplification_h
#define __vtkManifoldSimplification_h

#include "vtkSurface.h"
#include "RenderWindow.h"
#include "vtkPriorityQueue.h"

class VTK_EXPORT vtkManifoldSimplification : public vtkObject
{

public:

	// Create an instance of vtkManifoldSimplification
    static vtkManifoldSimplification *New();

	// Define the input mesh
	vtkSetObjectMacro(Input, vtkSurface)

	// Simplify the input mesh
	void Simplify();

	// Enable/Disable graphical display
	vtkSetMacro(Display, int)

	// Set the number of desired vertices
	vtkSetMacro(NumberOfOutputVertices, int)

protected:

	void AllocateMemory();
	void ReleaseMemory();
	void UpdateEdgePriority(vtkIdType Edge);

	RenderWindow *Window;

	// input mesh
	vtkSurface *Input;
	
	// the desired number of vertices
	int NumberOfOutputVertices;

	// priority queue used to compute geodesic voronoi regions
	vtkPriorityQueue *EdgesQueue;
	
	vtkBitArray *NonContractibleEdges;

	/// the constructor
	vtkManifoldSimplification();

	/// the desctructor
	~vtkManifoldSimplification();
	
	vtkIdList *FacesList;
	
	double *QuadricsArray;
	double **Quadrics;

	int Display;
};


#endif
