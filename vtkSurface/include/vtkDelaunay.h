/*=========================================================================

  Program:   Mesh to Delaunay Mesh
  Module:    vtkDelaunay.h
  Language:  C++
  Date:      2008/02
  Auteurs:   Sebastien VALETTE
             Arnaud Gelas

=========================================================================*/

/* ---------------------------------------------------------------------

* Copyright (c) CREATIS-LRMN (Centre de Recherche en Imagerie Medicale)
* Author : Sebastien Valette, Arnaud Gelas
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

#ifndef __vtkDelaunay_h
#define __vtkDelaunay_h

#include <vtkPriorityQueue.h>
#include "vtkSurface.h"
#include <vtkMath.h>

/**
 * This class converts a surface mesh into a mesh with delaunay properties by means of edge flips
 * See the paper "Delaunay Mesh Construction" from Dyer et al., SGP 2007
 */
class VTK_EXPORT vtkDelaunay : public vtkObject
{

public:

	// the constructor
	static vtkDelaunay *New();

	// Sets the mesh to modify
	void SetInputData(vtkSurface *Input);

	// This method flips the edges of the input mesh until all its edges meet the Delaunay criterion
	void DelaunayConform();

	// This method is intended to be used before vtkDelaunay::Update for each potentially flippable edge
	// this is faster than vtkDelaunay::DelaunayConform() which tests all the edges of the mesh
	void PushEdge(vtkIdType Edge);

	// After the user has pushed the edges that might be flipped, the call to Update() processes them
	void Update();
	
	// returns the number of flips performed during the call to Update()
	int GetNumberOfFlipedEdges()
	{return this->NumberOfPerformedFlips;};
	
	// Derive these methods when you want to do additionnal operations when an edge is flipped
	virtual void UpdateEdge(vtkIdType Edge){};
	virtual bool IsEdgeFlippable(vtkIdType Edge);

	static	double ComputeDelaunayConformingCriterion( vtkSurface *Mesh,const vtkIdType& Edge )
	{
		vtkIdType v1, v2, v3, v4;
		vtkIdType f1, f2;
		double p1[3], p2[3], p3[3], p4[3];
		double v31[3], v32[3], v41[3], v42[3];
		double dot3 = 0.;
		double dot4 = 0.;
		double norm31 = 0.;
		double norm32 = 0.;
		double norm41 = 0.;
		double norm42 = 0.;

		Mesh->GetEdgeFaces( Edge, f1, f2 );
		if (f2<0)
			return (-1.0);

		Mesh->GetEdgeVertices( Edge, v1, v2 );
		v3 = Mesh->GetThirdPoint( f1, v1, v2 );
		v4 = Mesh->GetThirdPoint( f2, v1, v2 );

		Mesh->GetPointCoordinates( v1, p1 );
		Mesh->GetPointCoordinates( v2, p2 );
		Mesh->GetPointCoordinates( v3, p3 );
		Mesh->GetPointCoordinates( v4, p4 );

		for( unsigned int dim = 0; dim < 3; dim++ )
		{
			v31[dim] = p1[dim] - p3[dim]; norm31 += v31[dim] * v31[dim];
			v32[dim] = p2[dim] - p3[dim]; norm32 += v32[dim] * v32[dim];
			dot3 += v31[dim] * v32[dim];

			v41[dim] = p1[dim] - p4[dim]; norm41 += v41[dim] * v41[dim];
			v42[dim] = p2[dim] - p4[dim]; norm42 += v42[dim] * v42[dim];
			dot4 += v41[dim] * v42[dim];
		}
		double den = norm31 * norm32;

		if( den != 0. )
		{
			den = sqrt( den );
			dot3 /= den;
		}

		if( dot3 > 1. )
			dot3 = 1.;
		if( dot3 < -1. )
			dot3 = -1.;

		den = norm41 * norm42;

		if( den != 0. )
		{
			den = sqrt( den );
			dot4 /= den;
		}

		if( dot4 > 1. )
			dot4 = 1.;
		if( dot4 < -1. )
			dot4 = -1.;

		return (acos( dot3 ) + acos( dot4 ) -  3.14159265358979323846);
	}

protected:

	vtkSurface *Mesh;

	vtkPriorityQueue *NonDelaunayEdges;
	
	int NumberOfPerformedFlips;

	vtkDelaunay();
	virtual	~vtkDelaunay();
};

#endif
