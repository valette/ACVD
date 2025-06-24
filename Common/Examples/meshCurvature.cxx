/*=========================================================================

Program:   meshCurvature : computes mesh hcurvature
Module:    vtkSurface
Language:  C++
Date:      2025/06
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

// .NAME meshCurvature
// .SECTION Description

#include "vtkSurface.h"
#include "vtkCurvatureMeasure.h"

int main( int argc, char *argv[] )
{

	if (argc<2) {
		cout<<"Usage : meshCurvature inputMesh outputMesh"<<endl;
		exit(1);
	}

	const char *outputFile = "output.vtp";
	vtkSurface *Mesh = vtkSurface::New();
	Mesh->CreateFromFile( argv[ 1 ] );
	if (argc > 2) {
		outputFile = argv[ 2 ];
	}
	cout << "Input mesh: " << argv[ 1 ] << endl;

	vtkCurvatureMeasure *curvatureMeasure = vtkCurvatureMeasure::New();
	curvatureMeasure->SetInputData( Mesh );
	curvatureMeasure->SetElementsType( 1 ); // 1 for points, 0 for cells
	curvatureMeasure->SetComputationMethod( 1 ); // 1 for polynomial fitting
	auto curvatureIndicator = curvatureMeasure->GetCurvatureIndicator()->GetItem( 0 );
	Mesh->GetPointData()->AddArray( curvatureIndicator );
	Mesh->GetPointData()->SetActiveScalars( curvatureIndicator->GetName() );
	Mesh->WriteToFile( outputFile );
	cout << "Output mesh: " << outputFile << endl;
}
