/*=========================================================================

Program:   mesh2off : a simple mesh conversion tool
Module:    vtkSurface
Language:  C++
Date:      2019/07
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

// .NAME mesh2vtk 
// .SECTION Description

#include "vtkSurface.h"
#include "vtkOFFWriter.h"

/// a simple mesh consersion tool
/// Usage : mesh2vtk inputfile
/// where inputfile is a mesh
/// the output is mesh.vtk

int main( int argc, char *argv[] )
{
	vtkSurface *Mesh;

	if (argc<2)
	{
		cout<<"Usage : mesh2off inputmesh"<<endl;
		exit(1);
	}

	// Load the mesh and create the vtkSurface data structure
	Mesh=vtkSurface::New();
	cout <<"load : "<<argv[1]<<endl;
	Mesh->CreateFromFile(argv[1]);
	vtkOFFWriter *Writer=vtkOFFWriter::New();
	Writer->SetInputData(Mesh);
	Writer->SetFileName("mesh.off");
	Writer->Write();
	cout<<"conversion to mesh.off finished!"<<endl;
	return (0);
}
