/*=========================================================================

Program:   meshDifference : computes coordinates differences on meshes with the same connectivity
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

// .NAME meshDifference 
// .SECTION Description

#include <vector>
#include "vtkSurface.h"
#include <vtkMath.h>
#include <vtkPolyDataWriter.h>


int main( int argc, char *argv[] )
{

	if (argc<2) {
		cout<<"Usage : meshDifference meshList.txt"<<endl;
		exit(1);
	}


	std::string line;
	std::ifstream file( argv[ 1 ] );
	std::vector<std::string> files;

	while( std::getline( file,line ) ) {

		files.push_back( line );
		cout << line << endl;

	}

	if ( files.size() % 2 ) {

		cout << "Eror : file list must contain an even number of files" << endl;
		exit ( 1 );

	}

	vtkDoubleArray *maxErrors = vtkDoubleArray::New();
	vtkDoubleArray *averageErrors = vtkDoubleArray::New();
	vtkDoubleArray *averageDisplacements = vtkDoubleArray::New();
	int numberOfPoints;
	double globalMaxDistance = 0;
	int globalMaxId = -1;

	for ( int i = 0; i < files.size(); i += 2 ) {

		vtkSurface *mesh1 = vtkSurface::New();
		vtkSurface *mesh2 = vtkSurface::New();
		mesh1->CreateFromFile( files[ i ].c_str() );
		mesh2->CreateFromFile( files[ i + 1 ].c_str() );

		if ( i == 0 ) {

			numberOfPoints = mesh1->GetNumberOfPoints();
			maxErrors->SetNumberOfValues( numberOfPoints );
			averageErrors->SetNumberOfValues( numberOfPoints );
			averageDisplacements->SetNumberOfComponents( 3 );
			averageDisplacements->SetNumberOfTuples( numberOfPoints );
			double d[ 3 ] = { 0, 0, 0};
			
			for ( int j = 0; j < numberOfPoints; j++ ) {
				maxErrors->SetValue( j, 0 );
				averageErrors->SetValue( j, 0 );
				averageDisplacements->SetTuple( j, d );
			}

		}

		if ( mesh1->GetNumberOfPoints() != numberOfPoints )
			cout << "error : " << files[ i ] << " should contain " << numberOfPoints << " vertices" << endl;

		if ( mesh1->GetNumberOfPoints() != numberOfPoints )
			cout << "error : " << files[ i + 1 ] << " should contain " << numberOfPoints << " vertices" << endl;


		double pt1[ 3 ], pt2[ 3 ], d[ 3 ];
		for ( int j = 0; j < numberOfPoints; j++ ) {

			mesh1->GetPointCoordinates( j, pt1 );
			mesh2->GetPointCoordinates( j, pt2 );
			double distance = sqrt( vtkMath::Distance2BetweenPoints( pt1, pt2 ) );
			if ( distance > maxErrors->GetValue( j ) ) maxErrors->SetValue( j, distance );
			averageErrors->SetValue( j, averageErrors->GetValue( j ) + distance / ( files.size() / 2 ) );
			averageDisplacements->GetTuple( j, d );
			for ( int k = 0; k < 3; k++ )
				d[ k ] += ( pt2[ k ] - pt1[ k ] )  / ( files.size() / 2 );
			averageDisplacements->SetTuple( j, d );

			if ( distance > globalMaxDistance ) {
				globalMaxDistance = distance;
				globalMaxId = i;
			}

		}

		mesh1->Delete();
		mesh2->Delete();

	}

	vtkSurface *mesh = vtkSurface::New();
	mesh->CreateFromFile( files[ 0 ].c_str() );
	mesh->GetPointData()->SetScalars( averageErrors );
	mesh->WriteToFile( "averageErrors.vtk" );
	mesh->Delete();

	mesh = vtkSurface::New();
	mesh->CreateFromFile( files[ 0 ].c_str() );
	mesh->GetPointData()->SetScalars( maxErrors );
	mesh->WriteToFile( "maxErrors.vtk" );
	mesh->Delete();

	mesh = vtkSurface::New();
	mesh->CreateFromFile( files[ 0 ].c_str() );
	mesh->GetPointData()->SetScalars( averageDisplacements );
	mesh->WriteToFile( "averageDisplacements.vtk" );
	mesh->Delete();

	cout << "Global Max distance = " << globalMaxDistance << " for id= " << globalMaxId << endl;

	return (0);
}
