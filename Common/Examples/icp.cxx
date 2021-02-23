/*=========================================================================

Program:   icp : iterative closest point
Module:    vtkSurface
Language:  C++
Date:      2021/02
Author:   Sebastien Valette

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

// .NAME icp 
// .SECTION Description

#include <vtkIterativeClosestPointTransform.h>
#include <vtkLandmarkTransform.h>
#include <vtkPLYWriter.h>
#include <vtkTransformPolyDataFilter.h>

#include "vtkSurface.h"

int main( int argc, char *argv[] ) {

	if  (argc < 3 ) {
		cout << "Usage : icp sourceMesh targetMesh" << endl;
		exit( 1 );
	}

	// Load the mesh and create the vtkSurface data structure
	vtkSurface *sourceMesh = vtkSurface::New();
	cout << "load source: " << argv[ 1 ] << endl;
	sourceMesh->CreateFromFile( argv[ 1 ] );

	vtkSurface *targetMesh = vtkSurface::New();
	cout << "load target: " << argv[ 2 ] << endl;
	targetMesh->CreateFromFile( argv[ 2 ] );

	vtkIterativeClosestPointTransform *icp = vtkIterativeClosestPointTransform::New();
	icp->SetMaximumNumberOfLandmarks( sourceMesh->GetNumberOfPoints() );
	cout << icp->GetNumberOfIterations () << " iterations" << endl;
	cout << icp->GetMaximumNumberOfLandmarks() << " landmarks maximum" << endl;
	icp->SetSource( sourceMesh );
	icp->SetTarget( targetMesh );

	icp->GetLandmarkTransform()->SetModeToRigidBody();
	icp->StartByMatchingCentroidsOn();
	icp->Update();
	cout << icp->GetNumberOfIterations () << " iterations" << endl;

	vtkTransformPolyDataFilter *transform = vtkTransformPolyDataFilter::New();
	transform->SetTransform( icp );
	transform->SetInputData( sourceMesh );
	transform->Update();

	vtkPLYWriter *Writer = vtkPLYWriter::New();
	Writer->SetInputData( transform->GetOutput() );
	Writer->SetFileName("output.ply");
	Writer->Write();

	//return (0);

}
