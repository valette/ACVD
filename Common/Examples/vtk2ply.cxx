/*=========================================================================

Program:   mesh2vtk : a simple mesh consersion tool
Module:    vtkSurface
Language:  C++
Date:      2011/07
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

#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPLYWriter.h>

/// a simple mesh consersion tool
/// Usage : vtk2ply inputfile
/// where inputfile is a mesh
/// the output is mesh.ply

int main( int argc, char *argv[] )
{
	if (argc<2)
	{
		cout<<"Usage : vtk2ply inputmesh"<<endl;
		exit(1);
	}

	cout <<"load : "<<argv[1]<<endl;

	vtkPolyDataReader *Reader=vtkPolyDataReader::New();
	Reader->SetFileName(argv[1]);
	Reader->Update();

	// Load the mesh and create the vtkSurface data structure
	vtkPLYWriter *Writer=vtkPLYWriter::New();
	Writer->SetInputData(Reader->GetOutput());
	Writer->SetFileName("mesh.ply");
	Writer->Write();
	cout<<"conversion to mesh.ply finished!"<<endl;
	return (0);
}
