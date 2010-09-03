/*=========================================================================

Program:   Tools for QEM processing
Module:    vtkSurface
Language:  C++
Date:      2007/10
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

#ifndef _VTKQUADRICTOOLS_H_
#define _VTKQUADRICTOOLS_H_

#define	DEFAULT_SV_THRESHOLD	0.001

#include <vtkMath.h>
#include <vtkObjectFactory.h>
#include <vtkTriangle.h>
#include "vtkSurface.h"

class VTK_EXPORT vtkQuadricTools : public vtkObject
{
public:

	//Create an instance of vtkQuadricTools
    static vtkQuadricTools *New();

	// Adds to Quadric the quadric equivalent to a plane passing by the given Point with the Given Normal
	// If FullQuadric=true then the whole quadric will be created (10 coefficients).
	// Otherwise, only the 9 first coefficients will be given
	static void AddPointWithNormalQuadric(double *Quadric, double *Point,
		double *Normal, double Factor=1.0, bool FullQuadric=true);

	// Adds the quadric computed from Face to Quadric, weighted by Factor
	// If FullQuadric=true then the whole quadric will be created (10 coefficients).
	// Otherwise, only the 9 first coefficients will be given
	static void AddTriangleQuadric(double *Quadric,vtkSurface *Mesh,int Face,bool FullQuadric=true);

	// Computes the displacement needed to reach the best position according to the quadric
	// MaxNumberOfUsedSingularValues defines the number of singular values used (generally 3)
	static int ComputeDisplacement(double *Quadric, double *Point, double *Displacement
				,int MaxNumberOfUsedSingularValues=3, double SVThreshold=DEFAULT_SV_THRESHOLD);

	// Projects the point on the position giving the minimum quadric error
	// returns the rank deficiency of the quadric.
	// MaxNumberOfUsedSingularValues defines the number of singular values used (generally 3)
	static int ComputeRepresentativePoint(double *Quadric, double *Point
					,int MaxNumberOfUsedSingularValues=3, double SVThreshold=DEFAULT_SV_THRESHOLD);

	static double Evaluate(double *Quadric, double *Point, bool FullQuadric=true);

	void GetPointQuadric(vtkSurface *Mesh, vtkIdType Vertex, double *Quadric, bool FullQuadric=true);

protected:

	// List used to get faces around a specific vertex
	vtkIdList *List;

	vtkQuadricTools();
	~vtkQuadricTools();
};

#endif
