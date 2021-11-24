/*=========================================================================

Program:   subdivideMesh : a simple mesh subdivision example
Module:    vtkSurface
Language:  C++
Date:      2021/11
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

// .NAME subdivideMesh 
// .SECTION Description

#include "vtkSurface.h"
#include "RenderWindow.h"

/// A simple example of mesh subdivision (mostly for benchmarking...)
/// Usage : subdivideMesh mesh n

int main( int argc, char *argv[] )
{

	// Load the mesh and create the vtkSurface data structure
	vtkSurface *Mesh = vtkSurface::New();
	cout <<"load : "<< argv[1]<<endl;
	Mesh->CreateFromFile(argv[1]);

	cout << Mesh->GetNumberOfPoints() << " Points" << endl;
	cout << Mesh->GetNumberOfCells() << " Faces" << endl;
	cout << Mesh->GetNumberOfEdges() << " Edges" << endl;

	for ( int i = 0; i < atoi( argv[ 2 ] ); i++ ) {
		cout << "*********************************" << endl;
		cout << "Subdividing..." << endl;
		Mesh = Mesh->Subdivide();
		cout << Mesh->GetNumberOfPoints() << " Points" << endl;
		cout << Mesh->GetNumberOfCells() << " Faces" << endl;
		cout << Mesh->GetNumberOfEdges() << " Edges" << endl;
	}

	return (0);
}
