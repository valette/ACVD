/*=========================================================================

Program:   Aproximated Centroidal Voronoi Diagrams (Quadric enhanced)
Module:    ACVD.cxx
Language:  C++
Date:      2006/04
Auteur:   Sebastien Valette,

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
// .NAME ACVDQ 
// .SECTION Description
#include <sstream>
#include <vtkPLYWriter.h>
#include <vtkSTLWriter.h>
#include <vtkCellData.h>

#include "vtkIsotropicDiscreteRemeshing.h"

/////////////////////////////////////////////////////////////////////////////////////////
// ACVDQ program:
/////////////////////////////////////////////////////////////////////////////////////////
// 
// Adaptive coarsening of triangular meshes (Quadrics enhanced)
// This program should be run with 3 arguments:
// run: "acvd file nvertices gradation [options]"
// file is the name of the mesh file to read
// nverticew is the desired number of vertices (note: if the number of input 
// gradation is the gradation parameter (0 is uniform, higher values give more and more importance 
//									to regions with high curvature)
//
// Additionnal options : 
// -d x : sets the graphics display (0 : no display. 1: display. 2 :iterative display)
//			default value : 1
//
// -s x : sets the subsampling threshold (Higher values give better results	 but the input 
//			mesh will be subdivided more times)
//			default value : 10
// -np x : sets the number of wanted threads (useful only with multi-processors machines)
//			default value : 1
//
// -o path : defines the output directory
//
//////////////////////////////////////////////////////////////////////////////////////////
// References:
// [1] " Approximated Centroidal Voronoi Diagrams for Uniform 
// Polygonal Mesh Coarsening", Valette & Chassery, Eurographics 2004.
// [2] "Adaptive Polygonal Mesh Simplification With Discrete Centroidal Voronoi Diagrams"
//  by, S. Valette, I. Kompatsiaris and J.-M. Chassery 
// + non published work....
/////////////////////////////////////////////////////////////////////////////////////////

int main( int argc, char *argv[] ) {

	//******************************************************************************************
	// Input parameters:
	int display = 0;				// defines whether there will be a graphic display (0: No, 1: yes)
	int numberOfSamples = 500;		// number of desired vertices
	double gradation = 0;			// gamma parameter for simplification (if gamma=0: uniform)
									// other appropriates values range between 0 and 2
	int subsamplingThreshold = 10;	// subsampling threshold
	char* outputDirectory = 0;		// output directory
	vtkIdList *fixedVertices = 0;
	//*******************************************************************************************

	char filename[ 5000 ];

	if( argc > 1 ) {

		cout << "load : " << argv[ 1 ] << endl;
		strcpy( filename, argv[ 1 ] );

	} else {

#ifdef __ACVDQP__
		cout<<"Usage : ACVDQP file nvertices gradation [options]"<<endl;
#else
		cout<<"Usage : ACVDQ file nvertices gradation [options]"<<endl;
#endif
		cout<<"nvertices is the desired number of vertices"<<endl;
		cout<<"gradation defines the influence of local curvature (0=uniform meshing)"<<endl;
		cout<<endl<<"Optionnal arguments : "<<endl;
		cout << "-b 0/1 : sets mesh boundary fixing off/on (default : 0)" << endl;
		cout<<"-s threshold : defines the subsampling threshold i.e. the input mesh will be subdivided until its number ";
		cout<<" of vertices is above nvertices*threshold (default=10)"<<endl;
		cout<<"-d 0/1/2 : enables display (default : 0)"<<endl;
		cout << "-l ratio : split the edges longer than ( averageLength * ratio )" << endl;
		cout << "-q 1/2/3 : qets number of eigenvalues used for quadric-based vertex relocation to 0/1/2 (default : 3)"<< endl;
		cout<<"-cd file : set custom imagedata file containing density information"<<endl;
		cout<<"-cmin value : set minimum custom indicator value"<<endl;
		cout<<"-cmax value : set maximum custom indicator value"<<endl;
		cout<<"-cf value : set custom indicator multiplication factor"<<endl;
		cout<<"-m 0/1 : enforce a manifold output ON/OFF (default : 0)"<<endl;
#ifdef DOmultithread
		cout<<"-np number : set the number of threads (default : number of cores)"<<endl;
#endif
		return 0;

	}

	vtkSurface *mesh = vtkSurface::New();
	vtkQIsotropicDiscreteRemeshing *remesh = vtkQIsotropicDiscreteRemeshing::New();
	mesh->CreateFromFile( filename );
	mesh->GetCellData()->Initialize();
	mesh->GetPointData()->Initialize();
	mesh->DisplayMeshProperties();

	// get mandatory arguments
	if( argc > 2 ) {

		numberOfSamples = atoi( argv[ 2 ] );

	} else {

		cout << "Number of vertices ? ";
		cin >> numberOfSamples;

	}

	if( argc > 3 ) gradation = atof( argv[ 3 ] );
	else {

		cout << "Gradation ? ";
		cin >> gradation;

	}

	// Parse optionnal arguments
	int argumentsIndex = 4;

	while ( argumentsIndex < argc ) {

		char* key = argv[ argumentsIndex ];
		char* value = argv[ argumentsIndex + 1 ];

		if ( strcmp( key, "-m" ) == 0 ) {

			remesh->SetForceManifold( atoi( value ) );
			cout << "Force Manifold=" << atoi( value ) << endl;

		} else if ( strcmp( key, "-s" ) == 0 ) {

			subsamplingThreshold = atoi( value );
			cout << "Subsampling Threshold=" << subsamplingThreshold << endl;

		} else if ( strcmp( key, "-d" ) == 0 ) {

			display = atoi( value );
			cout << "Display=" << display << endl;

		}

#ifdef DOmultithread
		if ( strcmp( key, "-np" ) == 0 )
		{
			int NumberOfThreads = atoi( value );
			cout << "Number of threads=" << NumberOfThreads << endl;
			remesh->SetNumberOfThreads( NumberOfThreads );
		}
#endif
		if ( strcmp( key, "-o" ) == 0 ) {

			outputDirectory = value;
			cout << "OutputDirectory: " << outputDirectory << endl;
			remesh->SetOutputDirectory( value );

		} else if ( strcmp( key, "-l" ) == 0 ) {

			cout << "Splitting edges longer than "
				<< atof( value ) << " times the average edge length" << endl;
			mesh->SplitLongEdges( atof( value ) );

		} else if ( strcmp( key,"-w" ) == 0 ) {

			cout << "Setting writing energy log file to " << atoi( value ) << endl;
			remesh->SetWriteToGlobalEnergyLog( atoi( value ) );

		}

#ifdef DOmultithread
		if ( strcmp( key, "-p" ) == 0 ) {

			cout << "Thread pooling ratio: " << atoi( value ) << endl;
			remesh->SetPoolingRatio( atoi( value ) );

		}
#endif

		if ( strcmp( key, "-q" ) == 0 ) {

			cout << "Setting number of eigenvalues for quadrics to " << atoi( value ) << endl;
			remesh->GetMetric()->SetQuadricsOptimizationLevel( atoi( value ) );

		} else if ( strcmp( key, "-cd" ) == 0 ) {

			cout << "Setting number custom file for density info : " << value << endl;
			remesh->SetInputDensityFile( value );

		} else if ( strcmp( key, "-cmax" ) == 0 ) {

			cout << "Setting maximum custom density to : " << value << endl;
			remesh->SetMaxCustomDensity( atof( value ) );

		} else if ( strcmp( key, "-cmin" ) == 0 ) {

			cout << "Setting minimum custom density to : " << value << endl;
			remesh->SetMinCustomDensity( atof( value ) );

		} if (strcmp( key, "-cf" ) == 0 ) {

			cout << "Setting custom density multiplication factor to : " << value << endl;
			remesh->SetCustomDensityMultiplicationFactor( atof( value ) );

		} else if ( strcmp( key, "-b" ) == 0 ) {

			cout << "Setting boundary fixing to : " << value << endl;
			remesh->SetBoundaryFixing( atoi( value ) );

		} else if ( strcmp( key, "-fv" ) == 0 ) {

			std::ifstream input;
			input.open( value );
			int id;
			fixedVertices = vtkIdList::New();
			while( input >> id ) fixedVertices->InsertNextId( id );
			input.close();

		} else if ( strcmp( key, "-ft" ) == 0 ) {

			std::ifstream input;
			input.open( value );
			bool* fixed = new bool[ mesh->GetNumberOfPoints() ];
			fixedVertices = vtkIdList::New();
			int id, n = 0;

			for ( int i = 0; i < mesh->GetNumberOfPoints(); i++ )
				fixed[ i ] = false;

			while( input >> id ) {

				n++;
				vtkIdType v1, v2, v3;
				mesh->GetFaceVertices( id, v1, v2, v3 );
				fixed[ v1 ] = fixed[ v2 ] = fixed[ v3 ] = true;

			}

			for ( int i = 0; i < mesh->GetNumberOfPoints(); i++ )
				if ( fixed[ i ] ) fixedVertices->InsertNextId( i );

			input.close();
			cout << "Added " << n << " constraints on triangles" << endl;
			delete [] fixed;

		}

		argumentsIndex += 2;

	}

	RenderWindow *window = 0;

	if ( display ) {

		window = RenderWindow::New();
		vtkPolyData *visu = vtkPolyData::New();
		visu->ShallowCopy( mesh );
		window->SetInputData( visu );
		visu->Delete();
		remesh->SetAnchorRenderWindow( window );
		window->Render();
		window->SetWindowName( filename );
		window->GetCamera()->Zoom( 1.2 );
		window->Interact();

	}

/*
	double bounds[ 6 ];
	mesh->GetBounds( bounds );
	double middle = 0.5 * ( bounds[ 0 ] + bounds [ 1 ] );
	cout << "middle : " << middle << endl;
	fixedVertices = vtkIdList::New();
	for ( int i = 0; i < mesh->GetNumberOfPoints(); i++ ) {
		double coords[ 3 ];
		mesh->GetPointCoordinates( i, coords );
		if ( coords[ 0 ] > middle ) fixedVertices->InsertNextId( i );
	}
	cout << "List size : " << fixedVertices->GetNumberOfIds() << endl;
*/

	remesh->SetInput( mesh );
	remesh->SetFileLoadSaveOption( 0 );
	remesh->SetConsoleOutput( 2 );
	remesh->SetSubsamplingThreshold( subsamplingThreshold );
	remesh->GetMetric()->SetGradation( gradation );
	remesh->SetDisplay( display );
	remesh->SetUnconstrainedInitialization( 1 );

	if ( fixedVertices ) {

		remesh->SetFixedClusters( fixedVertices );
		remesh->SetNumberOfClusters( numberOfSamples + fixedVertices->GetNumberOfIds() );
		cout << "Read " << fixedVertices->GetNumberOfIds() << " fixed Ids" << endl;

		for ( int i = 0; i < fixedVertices->GetNumberOfIds(); i++ )
			remesh->GetCluster( i )->AnchorItem = fixedVertices->GetId( i );

	} else remesh->SetNumberOfClusters( numberOfSamples );

	remesh->Remesh();

	// check that vertex constraints are respected
	if ( fixedVertices ) {

		vtkSurface *mesh2 = remesh->GetOutput();

		for ( int i = 0; i < fixedVertices->GetNumberOfIds(); i++ ) {

			double c1[ 3 ], c2[ 3 ];
			vtkIdType v = fixedVertices->GetId( i );
			mesh->GetPointCoordinates( v, c1 );
			mesh2->GetPointCoordinates( i, c2 );

			for ( int j = 0; j < 3; j++) {

				if ( c1[ j ] == c2[ j ] ) continue;
				cout << "Error, vertex " << v << " has been lost" << endl;
				exit( 1 );

			}


		}

		cout << "Constraints on vertices have been checked" << endl;
		fixedVertices->Delete();

	}

	// save the output mesh to .ply format
	char realFile[ 5000 ];
	
	if ( outputDirectory ) {

		strcpy ( realFile, outputDirectory );
		strcat ( realFile, "simplification.ply" );
	
	} else strcpy( realFile, "simplification.ply" );
	
	vtkPLYWriter *plyWriter = vtkPLYWriter::New();
	plyWriter->SetInputData( remesh->GetOutput() );
	plyWriter->SetFileName( realFile );
	plyWriter->Write();
	plyWriter->Delete();
	remesh->Delete();
	mesh->Delete();
	if ( display ) window->Delete();

}
