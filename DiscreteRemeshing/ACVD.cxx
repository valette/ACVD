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
#include <sstream>
#include <vtkPLYWriter.h>
#include <vtkSTLWriter.h>
#include <vtkCellData.h>

#include "vtkIsotropicDiscreteRemeshing.h"

/////////////////////////////////////////////////////////////////////////////////////////
// ACVD program:
/////////////////////////////////////////////////////////////////////////////////////////
// 
// Adaptive coarsening of triangular meshes
// This program should be run with 3 arguments:
// run: "acvd file nvertices gradation [options]"
// file is the name of the mesh file to read
// nvertices is the desired number of vertices (note: if the number of input 
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
/////////////////////////////////////////////////////////////////////////////////////////


int main( int argc, char *argv[] )
{

	//******************************************************************************************
	// Inside input parameters:
	int Display=1;				// defines whether there will be a graphic display (0: No, 1: yes)

	int NumberOfSamples=500;	// the number of desired vertices
	double Gradation=0;			// the gamma parameter for simplification (if gamma=0: uniform)
								// other appropriates values range between 0 and 2
	int SubsamplingThreshold=10;
	int QuadricsOptimizationLevel=1;
	char* OutputDirectory=0;		// the output directory
	//*******************************************************************************************

	char filename[500];

	if(argc>1)
	{
		cout <<"load : "<<argv[1]<<endl;
		strcpy(filename,argv[1]);
	}
	else
	{
		cout<<"Usage : ACVD file nvertices gradation [options]"<<endl;
		cout<<"nvertices is the desired number of vertices"<<endl;
		cout<<"gradation defines the influence of local curvature (0=uniform meshing)"<<endl;
		cout<<endl<<"Optionnal arguments : "<<endl;
		cout<<"-s threshold : defines the subsampling threshold i.e. the input mesh will be subdivided until its number ";
		cout<<" of vertices is above nvertices*threshold (default=10)"<<endl;
		cout<<"-d 0/1/2 : enables display (default : 0)"<<endl;
		cout<<"-q 0/1/2 : set the number of eigenvalues for quadrics post-processing (default : 3)"<<endl;
		cout<<"-cd file : set custom imagedata file containing density information"<<endl;
		cout<<"-cmin value : set minimum custom indicator value"<<endl;
		cout<<"-cmax value : set maximum custom indicator value"<<endl;
		cout<<"-cf value : set custom indicator multiplication factor"<<endl;
		cout<<"-m 0/1 : enforce a manifold output ON/OFF (default : 0)"<<endl;
		return (0);
	}

	vtkSurface *Mesh=vtkSurface::New();
	vtkIsotropicDiscreteRemeshing *Remesh=vtkIsotropicDiscreteRemeshing::New();

	Mesh->CreateFromFile(filename);
	Mesh->GetCellData()->Initialize();
	Mesh->GetPointData()->Initialize();

	Mesh->DisplayMeshProperties();

	// get mandatory arguments
	if(argc>2)
	{
		NumberOfSamples=atoi(argv[2]);
	}
	else
	{
		NumberOfSamples=3000;
		cout<<"Number of vertices ? ";
		cin>>NumberOfSamples;
	}

	if(argc>3)
	{
		Gradation=atof(argv[3]);
	}
	else
	{
		cout<<"Gradation ? ";
		cin>>Gradation;
	}

	// Parse optionnal arguments
	int ArgumentsIndex=4;
	while (ArgumentsIndex<argc)
	{
		if (strcmp(argv[ArgumentsIndex],"-m")==0)
		{
			Remesh->SetForceManifold(atoi(argv[ArgumentsIndex+1]));
			cout<<"Force Manifold="<<atoi(argv[ArgumentsIndex+1])<<endl;
		}

		if (strcmp(argv[ArgumentsIndex],"-s")==0)
		{
			SubsamplingThreshold=atoi(argv[ArgumentsIndex+1]);
			cout<<"Subsampling Threshold="<<SubsamplingThreshold<<endl;
		}

		if (strcmp(argv[ArgumentsIndex],"-d")==0)
		{
			Display=atoi(argv[ArgumentsIndex+1]);
			cout<<"Display="<<Display<<endl;
		}
#ifdef DOmultithread
		if (strcmp(argv[ArgumentsIndex],"-np")==0)
		{
			int NumberOfThreads=atoi(argv[ArgumentsIndex+1]);
			cout<<"Number of threads="<<NumberOfThreads<<endl;
			Remesh->SetNumberOfThreads(NumberOfThreads);
		}
#endif

		if (strcmp(argv[ArgumentsIndex],"-o")==0)
		{

			OutputDirectory=argv[ArgumentsIndex+1];
			cout<<"OutputDirectory: "<<OutputDirectory<<endl;
			Remesh->SetOutputDirectory(argv[ArgumentsIndex+1]);		
		}

		if (strcmp(argv[ArgumentsIndex],"-l")==0)
		{

			Mesh->SplitLongEdges(atof(argv[ArgumentsIndex+1]));
			cout<<"Splitting edges longer than "
			<<atof(argv[ArgumentsIndex+1])<<" times the average edge length"<<endl;
		}

		if (strcmp(argv[ArgumentsIndex],"-w")==0)
		{
			cout<<"Setting writing energy log file to "<<atoi(argv[ArgumentsIndex+1])<<endl;
			Remesh->SetWriteToGlobalEnergyLog(atoi(argv[ArgumentsIndex+1]));
		}

		if (strcmp(argv[ArgumentsIndex],"-q")==0)
		{
			cout<<"Setting number of eigenvalues for quadrics to "<<atoi(argv[ArgumentsIndex+1])<<endl;
			QuadricsOptimizationLevel=atoi(argv[ArgumentsIndex+1]);
		}

		if (strcmp(argv[ArgumentsIndex],"-cd")==0)
		{
			cout<<"Setting custom file for density info : "<<argv[ArgumentsIndex+1]<<endl;
			Remesh->SetInputDensityFile(argv[ArgumentsIndex+1]);
		}

		if (strcmp(argv[ArgumentsIndex],"-cmax")==0)
		{
			cout<<"Setting maximum custom density to : "<<argv[ArgumentsIndex+1]<<endl;
			Remesh->SetMaxCustomDensity(atof(argv[ArgumentsIndex+1]));
		}

		if (strcmp(argv[ArgumentsIndex],"-cmin")==0)
		{
			cout<<"Setting minimum custom density to : "<<argv[ArgumentsIndex+1]<<endl;
			Remesh->SetMinCustomDensity(atof(argv[ArgumentsIndex+1]));
		}

		if (strcmp(argv[ArgumentsIndex],"-cf")==0)
		{
			cout<<"Setting custom density multiplication factor to : "<<argv[ArgumentsIndex+1]<<endl;
			Remesh->SetCustomDensityMultiplicationFactor(atof(argv[ArgumentsIndex+1]));
		}

		ArgumentsIndex+=2;
	}

	RenderWindow *Window=0;
	if (Display!=0)
	{
		Window=RenderWindow::New();
		vtkPolyData *Visu=vtkPolyData::New();
		Visu->ShallowCopy(Mesh);
		Window->SetInput(Visu);
		Visu->Delete();
		Remesh->SetAnchorRenderWindow(Window);
		Window->Render();
		Window->SetWindowName(filename);
		Window->GetCamera()->Zoom(1.6);
		Window->Interact();
	}

	Remesh->SetInput(Mesh);
	Remesh->SetFileLoadSaveOption(0);
	Remesh->SetNumberOfClusters(NumberOfSamples);
	Remesh->SetConsoleOutput(2);
	Remesh->SetSubsamplingThreshold(SubsamplingThreshold);
	Remesh->GetMetric()->SetGradation(Gradation);

	if (Mesh->GetNumberOfPoints()>3000000)
		Remesh->SetDisplay(0);
	else
		Remesh->SetDisplay(Display);

	Remesh->Remesh();

	if (QuadricsOptimizationLevel!= 0)
	{
		// Note : this is an adaptation of Siggraph 2000 Paper : Out-of-core simplification of large polygonal models
		vtkIntArray *Clustering=Remesh->GetClustering();
		vtkPLYWriter *Writer=vtkPLYWriter::New();
		Writer->SetInput(Remesh->GetOutput());
		Writer->SetFileName("output_1.ply");
		Writer->Write();
		Writer->Delete();

		int Cluster,NumberOfMisclassedItems=0;

		double **ClustersQuadrics =new double*[NumberOfSamples];
		for (int i = 0; i < NumberOfSamples; i++)
		{
			ClustersQuadrics[i]=new double[9];
			for (int j=0;j<9;j++)
				ClustersQuadrics[i][j]=0;
		}

		vtkIdList *FList=vtkIdList::New();

		for (int i = 0; i < Remesh->GetNumberOfItems (); i++)
		{
			Cluster = Clustering->GetValue (i);
			if ((Cluster >= 0)&& (Cluster < NumberOfSamples))
			{
				if (Remesh->GetClusteringType() == 0)
				{
					vtkQuadricTools::AddTriangleQuadric(ClustersQuadrics[Cluster],Remesh->GetInput(),i,false);
				}
				else
				{
					Remesh->GetInput()->GetVertexNeighbourFaces(i,FList);
					for (int j=0;j<FList->GetNumberOfIds();j++)
						vtkQuadricTools::AddTriangleQuadric(ClustersQuadrics[Cluster]
								,Remesh->GetInput(),FList->GetId(j),false);				
				}
			}
			else
				NumberOfMisclassedItems++;
		}
		FList->Delete();

		if (NumberOfMisclassedItems)
			cout<<NumberOfMisclassedItems<<" Items with wrong cluster association"<<endl;

		double P[3];
		for (int i = 0; i < NumberOfSamples; i++)
		{
			Remesh->GetOutput()->GetPoint (i, P);
			vtkQuadricTools::ComputeRepresentativePoint(ClustersQuadrics[i], P,QuadricsOptimizationLevel);
			Remesh->GetOutput()->SetPointCoordinates (i, P);
			delete[] ClustersQuadrics[i];
		}
		delete [] ClustersQuadrics;

		Mesh->GetPoints()->Modified ();


		cout<<"After Quadrics Post-processing : "<<endl;
		Remesh->GetOutput()->DisplayMeshProperties();

		if (Display > 0)
		{
			RenderWindow *OptimizedMeshWindow=RenderWindow::New();
			OptimizedMeshWindow->AttachToRenderWindow(Remesh->GetDisplayWindow());
			OptimizedMeshWindow->SetInput(Remesh->GetOutput());
			OptimizedMeshWindow->SetWindowName("Coarsened model (quadric based placement)");
			OptimizedMeshWindow->Render ();
			OptimizedMeshWindow->Interact ();
		}
	}



	// save the output mesh to .ply format
	char REALFILE[500];
	if (OutputDirectory)
	{
		strcpy (REALFILE,OutputDirectory);
		strcat (REALFILE,"simplification.ply");
	}
	else
		strcpy(REALFILE,"simplification.ply");

	vtkPLYWriter *plyWriter=vtkPLYWriter::New();
	plyWriter->SetInput(Remesh->GetOutput());
	plyWriter->SetFileName(REALFILE);
	plyWriter->Write();
	plyWriter->Delete();
	Remesh->Delete();
	Mesh->Delete();
	if (Display!=0)
		Window->Delete();
}
