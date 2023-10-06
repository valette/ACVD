/*=========================================================================

Program:   Aproximated Centroidal Voronoi Diagrams
Module:    ACVD.cxx
Language:  C++
Date:      2003/11
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

// .NAME ACVD 
// .SECTION Description

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <iostream>

#include "vtkIsotropicDiscreteRemeshing.h"

using namespace std;

/////////////////////////////////////////////////////////////////////////////////////////
// ACVD program:
// 
// Adaptive coarsening of triangular meshes
// References:
// [1] " Approximated Centroidal Voronoi Diagrams for Uniform 
// Polygonal Mesh Coarsening", Valette & Chassery, Eurographics 2004.
// [2] "Adaptive Polygonal Mesh Simplification With Discrete Centroidal Voronoi Diagrams"
//  by, S. Valette, I. Kompatsiaris and J.-M. Chassery 
/////////////////////////////////////////////////////////////////////////////////////////


int main( int argc, char *argv[] )
{

	int Display=0;				// defines whether there will be a graphic display (0: No, 1: yes)

	int NumberOfSamples = 0;	// number of desired vertices
	double Gradation = 0;		// gamma parameter for simplification (if gamma=0: uniform)
								// other appropriates values range between 0 and 2
	int SubsamplingThreshold = 10;
	int QuadricsOptimizationLevel = 1;

	char* OutputDirectory = 0;
	char outputfile[500];

	strcpy (outputfile, "simplification.ply");

	if(argc > 3) {
		std::cout << "load : " << argv[1] << endl;
		NumberOfSamples = atoi(argv[2]);
		Gradation = atof(argv[3]);
	} else {

#ifdef __ACVDQP__
		cout << "Usage : ACVDP file nvertices gradation [options]" << endl;
#else
		cout << "Usage : ACVD file nvertices gradation [options]" << endl;
#endif
		cout << "nvertices is the desired number of vertices" << endl;
		cout << "gradation defines the influence of local curvature (0=uniform meshing)" << endl;
		cout << endl << "Optionnal arguments : " << endl;
		cout << "-b 0/1 : sets mesh boundary fixing off/on (default : 0)" << endl;
		cout << "-s threshold : defines the subsampling threshold i.e. the input mesh will be subdivided until inputVertices > outputVertices * ratio" << endl;
		cout << "-o directory : sets the output directory " << endl;
		cout << "-of file : sets the output file name " << endl;
		cout << " of vertices is above nvertices*threshold (default=10)" << endl;
		cout << "-d 0/1/2 : enables display (default : 0)" << endl;
		cout << "-l ratio : split the edges longer than ( averageLength * ratio )" << endl;
		cout << "-q 0/1/2 : set the number of eigenvalues for quadrics post-processing (default : 3)" << endl;
		cout << "-cd file : set custom imagedata file containing density information" << endl;
		cout << "-cmin value : set minimum custom indicator value" << endl;
		cout << "-cmax value : set maximum custom indicator value" << endl;
		cout << "-cf value : set custom indicator multiplication factor" << endl;
		cout << "-m 0/1 : enforce a manifold output ON/OFF (default : 0)" << endl;
#ifdef DOmultithread
		cout<<"-np number : set the number of threads (default : number of cores)"<<endl;
#endif
		return (0);
	}

	vtkSurface *Mesh = vtkSurface::New();
	Mesh->CreateFromFile(argv[1]);
	Mesh->GetCellData()->Initialize();
	Mesh->GetPointData()->Initialize();
	Mesh->DisplayMeshProperties();

	vtkIsotropicDiscreteRemeshing *Remesh = vtkIsotropicDiscreteRemeshing::New();

	// Parse optionnal arguments
	int ArgumentsIndex = 4;
	while (ArgumentsIndex < argc) {
		char* key = argv[ArgumentsIndex];
		char* value = argv[ArgumentsIndex + 1];

		if (strcmp(key, "-m") == 0) {
			Remesh->SetForceManifold(atoi(value));
			cout << "Force Manifold=" << atoi(value) << endl;
		}

		if (strcmp(key, "-s") == 0) {
			SubsamplingThreshold = atoi(value);
			cout << "Subsampling Threshold=" << SubsamplingThreshold << endl;
		}

		if (strcmp(key, "-d") == 0) {
			Display = atoi(value);
			cout << "Display=" << Display << endl;
		}
#ifdef DOmultithread
		if (strcmp(key, "-np") == 0) {
			int NumberOfThreads = atoi(value);
			cout << "Number of threads=" << NumberOfThreads << endl;
			Remesh->SetNumberOfThreads(NumberOfThreads);
		}
#endif

		if (strcmp(key, "-o") == 0) {
			OutputDirectory = value;
			cout << "OutputDirectory: " << OutputDirectory << endl;
			Remesh->SetOutputDirectory(value);		
		}

		if (strcmp(key, "-of") == 0) {
			strcpy(outputfile, value);
			cout << "Output file name: " << outputfile << endl;
		}

		if (strcmp(key, "-l") == 0) {
			cout << "Splitting edges longer than " <<
				atof(value) << " times the average edge length" << endl;
			Mesh->SplitLongEdges(atof(value));
		}

		if (strcmp(key, "-w") == 0) {
			cout << "Setting writing energy log file to " << atoi(value) << endl;
			Remesh->SetWriteToGlobalEnergyLog(atoi(value));
		}

		if (strcmp(key, "-q") == 0) {
			cout << "Setting number of eigenvalues for quadrics to " << atoi(value) << endl;
			QuadricsOptimizationLevel = atoi(value);
		}

		if (strcmp(key, "-cd") == 0) {
			cout << "Setting custom file for density info : " << value << endl;
			Remesh->SetInputDensityFile(value);
		}

		if (strcmp(key, "-cmax") == 0) {
			cout << "Setting maximum custom density to : " << value << endl;
			Remesh->SetMaxCustomDensity(atof(value));
		}

		if (strcmp(key,"-cmin")==0) {
			cout<<"Setting minimum custom density to : " << value << endl;
			Remesh->SetMinCustomDensity(atof(value));
		}

		if (strcmp(key, "-cf") == 0) {
			cout<<"Setting custom density multiplication factor to : "<< value << endl;
			Remesh->SetCustomDensityMultiplicationFactor(atof(value));
		}

		if (strcmp(key, "-b") == 0) {
			cout << "Setting boundary fixing to : " << value << endl;
			Remesh->SetBoundaryFixing(atoi(value));
		}
		ArgumentsIndex += 2;
	}

	RenderWindow *Window = 0;
	if (Display) {
		Window = RenderWindow::New();
		vtkPolyData *Visu=vtkPolyData::New();
		Visu->ShallowCopy(Mesh);
		Window->SetInputData(Visu);
		Visu->Delete();
		Remesh->SetAnchorRenderWindow(Window);
		Window->Render();
		Window->SetWindowName(argv[1]);
		Window->GetCamera()->Zoom(1.6);
		Window->Interact();
	}

	Remesh->SetInput(Mesh);
	Remesh->SetFileLoadSaveOption(0);
	Remesh->SetNumberOfClusters(NumberOfSamples);
	Remesh->SetConsoleOutput(2);
	Remesh->SetSubsamplingThreshold(SubsamplingThreshold);
	Remesh->GetMetric()->SetGradation(Gradation);
	Remesh->SetDisplay(Display);
	Remesh->Remesh();

	if (QuadricsOptimizationLevel != 0) {
		// Note : this is an adaptation of Siggraph 2000 Paper :
		// Out-of-core simplification of large polygonal models
		vtkIntArray *Clustering = Remesh->GetClustering();

		char REALFILE[5000];
		char FileBeforeProcessing[500];
		strcpy (FileBeforeProcessing,"smooth_");
		strcat (FileBeforeProcessing, outputfile);
		if (OutputDirectory) {
			strcpy (REALFILE, OutputDirectory);
			strcat (REALFILE, FileBeforeProcessing);
		} else {
			strcpy (REALFILE, FileBeforeProcessing);
		}

		Remesh->GetOutput()->WriteToFile(REALFILE);

		int Cluster,NumberOfMisclassedItems = 0;

		double **ClustersQuadrics = new double*[NumberOfSamples];
		for (int i = 0; i < NumberOfSamples; i++) {
			ClustersQuadrics[i] = new double[9];
			for (int j = 0; j < 9; j++) {
				ClustersQuadrics[i][j] = 0;
			}
		}

		vtkIdList *FList = vtkIdList::New();

		for (int i = 0; i < Remesh->GetNumberOfItems (); i++) {
			Cluster = Clustering->GetValue (i);
			if ((Cluster >= 0)&& (Cluster < NumberOfSamples)) {
				if (Remesh->GetClusteringType() == 0) {
					vtkQuadricTools::AddTriangleQuadric(
						ClustersQuadrics[Cluster], Remesh->GetInput(), i, false);
				} else {
					Remesh->GetInput()->GetVertexNeighbourFaces(i, FList);
					for (int j = 0;j < FList->GetNumberOfIds(); j++)
						vtkQuadricTools::AddTriangleQuadric(
							ClustersQuadrics[Cluster], Remesh->GetInput(), FList->GetId(j), false);
				}
			} else {
				NumberOfMisclassedItems++;
			}
		}
		FList->Delete();

		if (NumberOfMisclassedItems) {
			cout << NumberOfMisclassedItems << " Items with wrong cluster association" << endl;
		}

		double P[3];
		for (int i = 0; i < NumberOfSamples; i++) {
			Remesh->GetOutput()->GetPoint (i, P);
			vtkQuadricTools::ComputeRepresentativePoint(ClustersQuadrics[i], P, QuadricsOptimizationLevel);
			Remesh->GetOutput()->SetPointCoordinates (i, P);
			delete[] ClustersQuadrics[i];
		}
		delete [] ClustersQuadrics;

		Mesh->GetPoints()->Modified ();

		cout << "After Quadrics Post-processing : " << endl;
		Remesh->GetOutput()->DisplayMeshProperties();

		if (Display) {
			RenderWindow *OptimizedMeshWindow = RenderWindow::New();
			OptimizedMeshWindow->AttachToRenderWindow(Remesh->GetDisplayWindow());
			OptimizedMeshWindow->SetInputData(Remesh->GetOutput());
			OptimizedMeshWindow->SetWindowName("Coarsened model (quadric based placement)");
			OptimizedMeshWindow->Render ();
			OptimizedMeshWindow->Interact ();
		}
	}

	// save the output mesh to .ply format
	char REALFILE[500];
	strcpy (REALFILE, "");
	if (OutputDirectory) {
		strcpy (REALFILE, OutputDirectory);
		strcat (REALFILE, "/");
	}
	strcat (REALFILE, outputfile);

	Remesh->GetOutput()->WriteToFile(REALFILE);

	Remesh->Delete();
	Mesh->Delete();
	if (Display) {
	//	Window->Delete();
	}
}
