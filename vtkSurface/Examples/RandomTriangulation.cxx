/*=========================================================================

Program:   Random Triangulation generator
Module:    wavemesh.cxx
Language:  C++
Date:      2006/12
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

// .NAME meshviewer 
// .SECTION Description
#include <vtkPLYWriter.h>
#include "vtkSurface.h"
#include "RenderWindow.h"
#include "vtkRandomTriangulation.h"

int main( int argc, char *argv[] )
{
	if (argc<3)
	{
		cout<<"Random Triangulation generator"<<endl;
		cout<<"usage : \"RandomTriangulation numberofvertices type\" "<<endl;
		cout<<"available types :"<<endl;
		cout<<"0 : uniform plane "<<endl;
		cout<<"1 : non uniform plane"<<endl;
		cout<<"2 : plane made of 4 different regions with different densities"<<endl;
		cout<<"3 : half pipe"<<endl;
		cout<<"4 : half pipe cut on one corner"<<endl;
		return (1);
	}
	
	vtkRandomTriangulation *Triangulation=vtkRandomTriangulation::New();
	vtkSurface *Mesh=Triangulation->BuildRandomTriangulation(atoi(argv[1]),atoi(argv[2]));
	
	Mesh->DisplayMeshProperties();
	
	RenderWindow *Window=RenderWindow::New();
	Window->SetInputData(Mesh);
	Window->Render();
	Window->Interact();
	
	vtkPLYWriter *Writer=vtkPLYWriter::New();
	Writer->SetInputData(Mesh);
	Writer->SetFileName("mesh.ply");
	Writer->Write();
	Writer->Delete();
	
	Mesh->Delete();
	Window->Delete();	
	return (0);
}	
	
