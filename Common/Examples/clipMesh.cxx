/*=========================================================================

Program:   clipMesh : clip mesh using a plane
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

// .NAME clipMesh
// .SECTION Description

#include <vtkClipPolyData.h>
#include <vtkPlane.h>
#include <vtkPLYWriter.h>
#include <vtkMath.h>

#include "vtkSurface.h"
#include "RenderWindow.h"

int main( int argc, char *argv[] ) {

	if  (argc < 3 ) {
		cout << "Usage : clipMesh sourceMesh ox oy oz nx ny nz" << endl;
		exit( 1 );
	}

	int display = 1;
	int argumentsIndex = 2;

	while( argumentsIndex < argc) {

		char *key = argv[ argumentsIndex ];
		char *value = argv[ argumentsIndex + 1 ];

		if (strcmp(key, "-d") == 0) {
			display = atoi( value );
		}

		argumentsIndex+=2;

	}

	// Load the mesh and create the vtkSurface data structure
	vtkSurface *mesh = vtkSurface::New();
	cout << "load mesh: " << argv[ 1 ] << endl;
	mesh->CreateFromFile( argv[ 1 ] );

	double origin[ 3 ], normal[ 3 ];

	for ( int i = 0; i < 3; i++ ) {

		origin[ i ] = atof( argv[ 2 + i ] );
		normal[ i ] = atof( argv[ 5 + i ] );

	}

	float norm = vtkMath::Norm( normal );
	for ( int i = 0; i < 3; i++ ) normal[ i ] /= norm;

	vtkPlane *plane = vtkPlane::New();
	plane->SetOrigin( origin );
	plane->SetNormal( normal );

	vtkClipPolyData *clip = vtkClipPolyData::New();
	clip->SetClipFunction( plane );
	clip->SetInputData( mesh );
	clip->Update();

	vtkPLYWriter *Writer = vtkPLYWriter::New();
	Writer->SetInputData( clip->GetOutput() );
	Writer->SetFileName( "output.ply" );
	Writer->Write();

	if ( display == 0 ) return( 0 );

	RenderWindow *window = RenderWindow::New();
	window->SetInputData( clip->GetOutput() );
	window->Interact();

}
