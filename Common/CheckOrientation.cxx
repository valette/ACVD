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

#include <vtkCommand.h>
#include <vtkCellData.h>
#include <vtkPointData.h>

#include "vtkSurface.h"
#include <vtkPolyDataNormals.h>
#include <vtkPolyDataWriter.h>

int main( int argc, char *argv[] )
{
	if (argc < 2) {
		cout << "Usage : CheckOrientation input" << endl;
		exit (1);
	}
	vtkSurface *Mesh = vtkSurface::New();
	Mesh->CreateFromFile( argv[1] );

	vtkPolyDataNormals *Normals = vtkPolyDataNormals::New();
	Normals->SetInputData (Mesh);
	Normals->SplittingOff();
	Normals->FlipNormalsOn();
	Normals->Update();

	vtkPolyData *Output = Normals->GetOutput();
	Output->GetPointData()->SetScalars(0);
	Output->GetCellData()->SetScalars(0);

	vtkPolyDataWriter *Writer = vtkPolyDataWriter::New();
	Writer->SetInputData(Output);
	Writer->SetFileName("good_orientation.vtk");
	Writer->Write();
	Writer->Delete();
	Normals->Delete();
	Mesh->Delete();
	return(0);
}
