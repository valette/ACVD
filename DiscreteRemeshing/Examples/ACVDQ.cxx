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


int main( int argc, char *argv[] )
{

	//******************************************************************************************
	// Inside input parameters:
	int Display=1;				// defines whether there will be a graphic display (0: No, 1: yes)


	int NumberOfSamples=500;	// the number of desired vertices
	double Gradation=0;			// the gamma parameter for simplification (if gamma=0: uniform)
								// other appropriates values range between 0 and 2
	int SubsamplingThreshold=10;
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
		cout << "-b 0/1 : sets mesh boundary fixing off/on (default : 0)" << endl;
		cout<<"-s threshold : defines the subsampling threshold i.e. the input mesh will be subdivided until its number ";
		cout<<" of vertices is above nvertices*threshold (default=10)"<<endl;
		cout<<"-d 0/1/2 : enables display (default : 0)"<<endl;
		cout<<"-cd file : set custom imagedata file containing density information"<<endl;
		cout<<"-cmin value : set minimum custom indicator value"<<endl;
		cout<<"-cmax value : set maximum custom indicator value"<<endl;
		cout<<"-cf value : set custom indicator multiplication factor"<<endl;
		cout<<"-m 0/1 : enforce a manifold output ON/OFF (default : 0)"<<endl;
		cout<<"-sf spare_factor : sets the spare factor"<<endl;
		return (0);
	}

	vtkSurface *Mesh=vtkSurface::New();
	vtkQIsotropicDiscreteRemeshing *Remesh=vtkQIsotropicDiscreteRemeshing::New();

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
		char* key = argv[ArgumentsIndex];
		char* value = argv[ArgumentsIndex + 1];

		if (strcmp(key,"-m")==0)
		{
			Remesh->SetForceManifold(atoi(value));
			cout<<"Force Manifold="<<atoi(value)<<endl;
		}

		if (strcmp(key,"-s")==0)
		{
			SubsamplingThreshold=atoi(value);
			cout<<"Subsampling Threshold="<<SubsamplingThreshold<<endl;
		}

		if (strcmp(key,"-d")==0)
		{
			Display=atoi(value);
			cout<<"Display="<<Display<<endl;
		}

#ifdef DOmultithread
		if (strcmp(key,"-np")==0)
		{
			int NumberOfThreads=atoi(value);
			cout<<"Number of threads="<<NumberOfThreads<<endl;
			Remesh->SetNumberOfThreads(NumberOfThreads);
		}
#endif
		if (strcmp(key,"-o")==0)
		{

			OutputDirectory=value;
			cout<<"OutputDirectory: "<<OutputDirectory<<endl;
			Remesh->SetOutputDirectory(value);
		}

		if (strcmp(key,"-l")==0)
		{

			Mesh->SplitLongEdges(atof(value));
			cout<<"Splitting edges longer than "
			<<atof(value)<<" times the average edge length"<<endl;
		}
		if (strcmp(key,"-w")==0)
		{
			cout<<"Setting writing energy log file to "<<atoi(value)<<endl;
			Remesh->SetWriteToGlobalEnergyLog(atoi(value));
		}
#ifdef DOmultithread
		if (strcmp(key,"-p")==0)
		{
			cout<<"Thread pooling ratio: "<<atoi(value)<<endl;
			Remesh->SetPoolingRatio(atoi(value));
		}
#endif

		if (strcmp(key,"-q")==0)
		{
			cout<<"Setting number of eigenvalues for quadrics to "<<atoi(value)<<endl;
			Remesh->GetMetric()->SetQuadricsOptimizationLevel(atoi(value));
		}

		if (strcmp(key,"-cd")==0)
		{
			cout<<"Setting number custom file for density info : "<<value<<endl;
			Remesh->SetInputDensityFile(value);
		}

		if (strcmp(key,"-cmax")==0)
		{
			cout<<"Setting maximum custom density to : "<<value<<endl;
			Remesh->SetMaxCustomDensity(atof(value));
		}

		if (strcmp(key,"-cmin")==0)
		{
			cout<<"Setting minimum custom density to : "<<value<<endl;
			Remesh->SetMinCustomDensity(atof(value));
		}

		if (strcmp(key,"-cf")==0)
		{
			cout<<"Setting custom density multiplication factor to : "<<value<<endl;
			Remesh->SetCustomDensityMultiplicationFactor(atof(value));
		}

		if (strcmp(key,"-sc")==0)
		{
			cout<<"Setting number of spare clusters to : "<<value<<endl;
			Remesh->SetMinNumberOfSpareClusters(atoi(value));
		}

		if (strcmp(key,"-sf")==0)
		{
			cout<<"Setting spare factor to : "<<value<<endl;
			Remesh->SetSpareFactor(atof(value));
		}

		if (strcmp(key, "-b") == 0) {
			cout << "Setting boundary fixing to : " << value << endl;
			Remesh->SetBoundaryFixing(atoi(value));
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
		Window->GetCamera()->Zoom(1.2);
		Window->Interact();
	}

	Remesh->SetInput(Mesh);
	Remesh->SetFileLoadSaveOption(0);
	Remesh->SetNumberOfClusters(NumberOfSamples);
	Remesh->SetConsoleOutput(2);
	Remesh->SetSubsamplingThreshold(SubsamplingThreshold);
	Remesh->GetMetric()->SetGradation(Gradation);
	Remesh->SetDisplay(Display);
	Remesh->SetConstrainedInitialization(1);
	Remesh->Remesh();

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
	plyWriter->SetInputData(Remesh->GetOutput());
	plyWriter->SetFileName(REALFILE);
	plyWriter->Write();
	plyWriter->Delete();
	Remesh->Delete();
	Mesh->Delete();
	if (Display!=0)
		Window->Delete();	
}
