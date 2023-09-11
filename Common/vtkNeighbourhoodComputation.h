/***************************************************************************
vtkNeighbourhoodComputation.h  -  description
-------------------
begin                : January 4 2006
copyright            : (C) 2006 by Sebastien Valette
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

#ifndef _VTKNEIGHBOURHOODCOMPUTATION_H_
#define _VTKNEIGHBOURHOODCOMPUTATION_H_

#include "vtkCommand.h"
#include "vtkObjectFactory.h"
#include "vtkSurface.h"
#include "vtkTag.h"

/// A Class to compute Neighbourhoods on meshes
/// The class is designed for efficient and possibly multithreaded computation
class VTK_EXPORT vtkNeighbourhoodComputation : public vtkObject
{
public:
	/// the public constructor
	static vtkNeighbourhoodComputation* New();

	/// Method to initialize the class for a specific Input
	void SetInputData(vtkSurface *Mesh);

	/// Compute the NRing around the Cell.
	void ComputeNRingCells(vtkIdType Cell,int RingSize,vtkIdList *FList);

	/// Compute the Cells around the input cell within the input Distance.
	void ComputeDistanceRingCells (vtkIdType Cell,double Distance, vtkIdList *FList);

	// Defines the type of the origin cells (0=faces 1=Vertices)
	void SetCellType (int Type)
	{ this->CellType=Type;};

	/// Returns the Input mesh
	vtkSurface *GetInput(){return (this->Input);};

protected:

	vtkNeighbourhoodComputation();
	~vtkNeighbourhoodComputation();

private:

	// Type of the origin cells (0=faces 1=Vertices)
	int CellType;
	
	// The input mesh
	vtkSurface *Input;

	// those arrays define which elements have already been visited
	vtkTag *VisitedCells;
	vtkTag *VisitedVertices;

	// IdLists statically created to speed up the neighborhood computation
	vtkIdList *VList;
	vtkIdList *VList2;
	vtkIdList *EList;
};
#endif
