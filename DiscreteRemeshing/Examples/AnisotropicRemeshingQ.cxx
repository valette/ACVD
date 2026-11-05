/*=========================================================================

Program:   Aproximated Centroidal Voronoi Diagrams
Module:    ACVD.cxx
Language:  C++
Date:      2003/11
Auteur:   Sebastien Valette,

=========================================================================*/
// .NAME AnisotropicRemeshingQ 
// .SECTION Description
#include <sstream>
#include <vtkSMPTools.h>
#include <vtkPLYWriter.h>
#include <vtkSTLWriter.h>
#include <vtkCellData.h>

#include "vtkDiscreteRemeshing.h"

#include "vtkIsotropicMetricForClustering.h"
#include "vtkQuadricAnisotropicMetricForClustering.h"
#include "vtkQEMetricForClustering.h"



#include "vtkTrianglesProcessing.h"
#include "vtkVerticesProcessing.h"

// the remesher passes itself to vtkDiscreteRemeshing as its most-derived class
template < class Metric > class VerticesRemeshing :
	public vtkVerticesProcessing < vtkDiscreteRemeshing < Metric, VerticesRemeshing < Metric > > >
{
public:
	static VerticesRemeshing *New() { return new VerticesRemeshing; }
};

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
	// VTK builds default to the Sequential SMP backend; STDThread is always available
	vtkSMPTools::SetBackend("STDThread");

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

	typedef VerticesRemeshing<vtkIsotropicMetricForClustering> IsotropicRemeshing;

	typedef VerticesRemeshing<vtkQEMetricForClustering> QEMRemeshing;

//	typedef VerticesRemeshing<vtkL21MetricForClustering> L21Remeshing;

	typedef VerticesRemeshing<vtkQuadricAnisotropicMetricForClustering> QuadricAnisotropicRemeshing;


	QuadricAnisotropicRemeshing *Remesh=QuadricAnisotropicRemeshing::New();

	if(argc>1)
	{
		cout <<"load : "<<argv[1]<<endl;
		strcpy(filename,argv[1]);
	}
	else
	{
#ifdef __ANISOTROPICREMESHINGQP__
		cout<<"Usage : AnisotropicRemeshingQ file nvertices gradation [options]"<<endl;
#else
		cout<<"Usage : AnisotropicRemeshingQ file nvertices gradation [options]"<<endl;
#endif
		cout<<"nvertices is the desired number of vertices"<<endl;
		cout<<"gradation defines the influence of local curvature (0=uniform meshing)"<<endl;
		cout<<endl<<"Optionnal arguments : "<<endl;
		cout << "-b 0/1 : sets mesh boundary fixing off/on (default : 0)" << endl;
		cout<<"-d 0/1/2 : enables display (default : 0)"<<endl;
		cout << "-l ratio : split the edges longer than ( averageLength * ratio )" << endl;
		cout << "-b 0/1 : sets mesh boundary fixing off/on (default : 0)" << endl;
		cout << "-q 1/2/3 : qets number of eigenvalues used for quadric-based vertex relocation to 0/1/2 (default : 3)"<< endl;
#ifdef DOmultithread
		cout<<"-np number : set the number of threads (default : number of cores)"<<endl;
#endif
		return (0);
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
			vtkSMPTools::Initialize(NumberOfThreads);
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

		if (strcmp(argv[ArgumentsIndex],"-q")==0)
		{
			cout<<"Setting number of eigenvalues for quadrics to "<<atoi(argv[ArgumentsIndex+1])<<endl;
			Remesh->GetMetric()->SetQuadricsOptimizationLevel(atoi(argv[ArgumentsIndex+1]));
		}

		if (strcmp(argv[ArgumentsIndex], "-b") == 0) {
			cout << "Setting boundary fixing to : " << argv[ArgumentsIndex+1] << endl;
			Remesh->SetBoundaryFixing(atoi(argv[ArgumentsIndex+1]));
		}

		ArgumentsIndex+=2;
	}

	RenderWindow *Window;
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

	Remesh->SetInput(Mesh);
	Remesh->SetNumberOfClusters(NumberOfSamples);
	Remesh->SetConsoleOutput(2);
	Remesh->SetSubsamplingThreshold(SubsamplingThreshold);
	Remesh->GetMetric()->SetGradation(Gradation);

//	Remesh->SetInitialSamplingType(0);
//	Remesh->SetDisplayCharacteristicsWindowOn();
//	Remesh->SetAlgorithmType(1);

	Remesh->SetDisplay(Display);
	Remesh->Remesh();
	
	// save the output mesh to .ply format
	char REALFILE[500];
	if (OutputDirectory)
	{
		strcpy (REALFILE,OutputDirectory);
		strcat (REALFILE,"Remeshing.ply");
		cout<<"OutputDirectory: "<<OutputDirectory<<endl;
	}
	else
		strcpy(REALFILE,"Remeshing.ply");

	vtkPLYWriter *plyWriter=vtkPLYWriter::New();
	plyWriter->SetInputData(Remesh->GetOutput());
	plyWriter->SetFileName(REALFILE);
	plyWriter->Write();
	plyWriter->Delete();
}
