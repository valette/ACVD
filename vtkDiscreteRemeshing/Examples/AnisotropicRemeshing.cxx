/*=========================================================================

Program:   Aproximated Centroidal Voronoi Diagrams
Module:    ACVD.cxx
Language:  C++
Date:      2003/11
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

// .NAME AnisotropicRemeshing 
// .SECTION Description
#include <sstream>
#include <vtkPLYWriter.h>
#include <vtkSTLWriter.h>
#include <vtkCellData.h>

#include "vtkDiscreteRemeshing.h"

#include "vtkIsotropicMetricForClustering.h"
#include "vtkAnisotropicMetricForClustering.h"
#include "vtkQEMetricForClustering.h"
#include "vtkL21MetricForClustering.h"


#include "vtkTrianglesProcessing.h"
#include "vtkVerticesProcessing.h"

//#include "vtkSubdivisionRemeshing.h"

/////////////////////////////////////////////////////////////////////////////////////////
// 
// Adaptive coarsening of triangular meshes
// This program should be run with 3 arguments:
// run: "acvd file nvertices gradation [options]"
// arg1 is the name of the mesh file to read
// arg2 is the desired number of vertices (note: if the number of input 
// arg3 is the gradation parameter (0 is uniform, higher values give more and more importance 
//									to regions with high curvature)
//
// Additionnal options : 
// -d x : sets the graphics display (0 : no display. 1: display. 2 :iterative display)
//			default value : 0
//
// -s x : sets the subsampling threshold (Higher values give better results	 but the input 
//			mesh will be subdivided more times)
//			default value : 10
// -np x : sets the number of wanted processes (useful only with multi-processors machines)
//
//////////////////////////////////////////////////////////////////////////////////////////

int main( int argc, char *argv[] )
{

//******************************************************************************************
	// Inside input parameters:
	int Display=0;			// defines whether there will be a graphic
					// display (0: No, 1: yes)

	int NumberOfSamples=200;	// the number of desired vertices

	double Gradation=0;		// the gamma parameter for simplification 
					// (if gamma=0: uniform)
					// other appropriates values range between 0 and 2

	int SubsamplingThreshold=10;
	char* OutputDirectory=0;		// the output directory 

//*******************************************************************************************

	char filename[500];

	vtkSurface *Mesh=vtkSurface::New();

	typedef vtkDiscreteRemeshing<vtkIsotropicMetricForClustering> IsotropicRemeshing;

	typedef vtkDiscreteRemeshing<vtkQEMetricForClustering> QEMRemeshing;

	typedef vtkDiscreteRemeshing<vtkL21MetricForClustering> L21Remeshing;

	typedef vtkDiscreteRemeshing<vtkAnisotropicMetricForClustering> AnisotropicRemeshing;

	vtkVerticesProcessing<AnisotropicRemeshing> *Remesh=
		vtkVerticesProcessing<AnisotropicRemeshing>::New();

	if(argc>1)
	{
		cout <<"load : "<<argv[1]<<endl;
		strcpy(filename,argv[1]);
	}
	else
	{
		strcpy(filename,"../Datas/hand.vtk");
		strcpy(filename,"/home/sebastien/CVSROOT/datas/horse.vtk");
	}

	Mesh->CreateFromFile(filename);
	Mesh->GetCellData()->Initialize();
	Mesh->GetPointData()->Initialize();

	Mesh->DisplayMeshProperties();

	if (0)
	{
		vtkPLYWriter *plyWriter=vtkPLYWriter::New();
		plyWriter->SetInputData(Mesh);
		plyWriter->SetFileName("input.ply");
		plyWriter->Write();
		plyWriter->Delete();

	}

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
	cout<<argc<<" Arguments"<<endl;
	while (ArgumentsIndex<argc)
	{
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
		}

		if (strcmp(argv[ArgumentsIndex],"-l")==0)
		{

			Mesh->SplitLongEdges(atof(argv[ArgumentsIndex+1]));
			cout<<"Splitting edges longer than "
			<<atof(argv[ArgumentsIndex+1])<<" times the average edge length"<<endl;
		}

		ArgumentsIndex+=2;
	}


	// Display the mesh

	RenderWindow *Window=0;
	if (Display!=0)
	{
		Window=RenderWindow::New();
		vtkPolyData *Visu=vtkPolyData::New();
		Visu->ShallowCopy(Mesh);
		Window->SetInputData(Visu);
		Remesh->SetAnchorRenderWindow(Window);
		Window->Render();
		Window->SetWindowName(filename);
		Window->GetCamera()->Zoom(1.6);
		Window->Interact();
	}


	// Remesh

	Remesh->SetInput(Mesh);
	Remesh->SetNumberOfClusters(NumberOfSamples);
	Remesh->SetConsoleOutput(2);
	Remesh->SetSubsamplingThreshold(SubsamplingThreshold);
	Remesh->GetMetric()->SetGradation(Gradation);

	if (Mesh->GetNumberOfPoints()>3000000)
		Remesh->SetDisplay(0);
	else
		Remesh->SetDisplay(Display);

	Remesh->Remesh();

	// save the output mesh to .ply format
	char REALFILE[500];
	if (OutputDirectory)
	{
		strcpy (REALFILE,OutputDirectory);
		strcat (REALFILE,"output.ply");
		cout<<"OutputDirectory: "<<OutputDirectory<<endl;
	}
	else
		strcpy(REALFILE,"output.ply");

	vtkPLYWriter *plyWriter=vtkPLYWriter::New();
	plyWriter->SetInputData(Remesh->GetOutput());
	plyWriter->SetFileName(REALFILE);
	plyWriter->Write();
	plyWriter->Delete();
	if (Display!=0)
		Window->Delete();
		
	Mesh->Delete();
	Remesh->Delete();	
}
