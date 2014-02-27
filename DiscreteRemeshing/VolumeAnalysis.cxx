/*=========================================================================

Program:   VolumeAnalysis : generate meshes according to volume labels
Module:    vtkSurface
Language:  C++
Date:      2011/02
Auteur:   Sebastien Valette

=========================================================================*/
// .NAME VolumeAnalysis 
// .SECTION Description

#include <sstream>
#include <vector>
#include <math.h>
#include <map>

#include <vtkImageData.h>
#include <vtkMetaImageReader.h>
#include <vtkIdList.h>
#include <vtkPolyDataWriter.h>
#include <vtkPLYWriter.h>
#include <vtkSTLWriter.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkMultiThreader.h>
#include <vtkXMLDataElement.h>
#include <vtkXMLUtilities.h>
#include <vtkTimerLog.h>
#include <vtkMutexLock.h>
#include <vtkQuadricClustering.h>
#include <vtkImageResample.h>
#include <vtkDiscreteMarchingCubes.h>
#include <vtkMINCImageReader.h>
#include <vtkWindowedSincPolyDataFilter.h>
#include <vtkSortDataArray.h>
#include <vtkFillHolesFilter.h>
#include <vtkImageCast.h>

#include "vtkIsotropicDiscreteRemeshing.h"
#include "vtkAnisotropicDiscreteRemeshing.h"

#define NumberOfTimingTypes 6

class MyThreaderHelperClass
{
public:
    vtkImageData *Image;
    vtkIdList *Labels;
    vtkMutexLock *Lock;
    int MaximumNumberOfVertices;
	std::vector < std::vector<double> > ProcessingTimes;
	std::vector < std::vector<double> > MaximumTimes;
	int SimplificationType;
	int NumberOfSmoothingSteps;
	char OutputDirectory[5000];
	double Gradation;
	char OutputFormat[10];
	int ForceManifold;
	int FillHoles;
	int Anisotropy;
	int KeepBiggestComponent;
	std::map<int, std::string> NamesMap;

	void SetOutputFormat(char *Format)
	{
		char LowerCaseFormat[10];
		strcpy(LowerCaseFormat,Format);

		if (LowerCaseFormat != NULL)
		{
			char *p;
			for (p = LowerCaseFormat; *p; ++p)
				*p = tolower(*p);

			strcpy(this->OutputFormat, LowerCaseFormat);
		}
	}

	char *GetOutputFormat()
	{
		return (this->OutputFormat);
	}

	MyThreaderHelperClass()
	{
		this->KeepBiggestComponent=0;
		this->ForceManifold=0;
		this->FillHoles=0;
		this->Anisotropy=0;
		this->Lock=vtkMutexLock::New();
		this->MaximumNumberOfVertices=(1<<15) -1;
		this->SimplificationType=0;
		this->NumberOfSmoothingSteps=0;
		strcpy (this->OutputDirectory,"");
		this->Gradation=0;
		this->SetOutputFormat((char *)"vtk");
	}

	~MyThreaderHelperClass()
	{
		this->Lock->Delete();
	}
};

VTK_THREAD_RETURN_TYPE ThreadedSurfaceExtraction (void *arg)
{
	vtkMultiThreader::ThreadInfo *Info = (vtkMultiThreader::ThreadInfo*) arg;
	MyThreaderHelperClass *Helper = (MyThreaderHelperClass *) Info->UserData;
	int MyId = Info->ThreadID;
	int NumberOfThreads=Info->NumberOfThreads;
	vtkIdList *Labels=Helper->Labels;

	vtkTimerLog *Timer=vtkTimerLog::New();

	Helper->Lock->Lock();
	if (Helper->ProcessingTimes.size()==0)
	{
		Helper->ProcessingTimes.resize(NumberOfThreads);
		Helper->MaximumTimes.resize(NumberOfThreads);
	}
	for (int i=0;i<NumberOfTimingTypes;i++)
	{
		Helper->ProcessingTimes[MyId].push_back(0);
		Helper->MaximumTimes[MyId].push_back(0);
	}
	Helper->Lock->Unlock();

	vtkImageData *Image=vtkImageData::New();
	Image->ShallowCopy(Helper->Image);

	// create the polydataWriter in correct format
	vtkPolyDataWriter *Writer=0;
	char *Format=Helper->GetOutputFormat();
	char extension[4];
	strcpy (extension,"vtk");
	if (strstr(Format,extension)!=NULL)
		Writer=vtkPolyDataWriter::New();
	else
	{
		strcpy (extension,"ply");
		if (strstr(Format,extension)!=NULL)
			Writer=(vtkPolyDataWriter*) vtkPLYWriter::New();
		else
		{
			strcpy (extension,"stl");
			if (strstr(Format,extension)!=NULL)
				Writer=(vtkPolyDataWriter*) vtkSTLWriter::New();		
		}
	}

	vtkDiscreteMarchingCubes *Contour=vtkDiscreteMarchingCubes::New();
	Contour->ComputeNormalsOff();
	Contour->ComputeGradientsOff ();
	Contour->ComputeScalarsOff();
	Contour->SetInputData(Image);
	for (int i=MyId; i<Labels->GetNumberOfIds(); i+=NumberOfThreads)
	{
		double StartTime=Timer->GetUniversalTime();
		vtkIdType Level=Labels->GetId(i);
		Timer->StartTimer();
		Contour->SetValue(0,Level);
		Contour->Update();
		vtkPolyData *Mesh=Contour->GetOutput();

		if (Helper->KeepBiggestComponent!=0)
		{
			vtkSurface *ToClean=vtkSurface::New();
			ToClean->CreateFromPolyData(Mesh);
			vtkSurface *Cleaned=ToClean->GetBiggestConnectedComponent();
			Cleaned->EnsureOutwardsNormals();
			Mesh->ShallowCopy(Cleaned);
			Cleaned->Delete();
			ToClean->Delete();
			/*
			vtkPolyDataConnectivityFilter *Connectivity=vtkPolyDataConnectivityFilter::New();
			Connectivity->SetExtractionModeToLargestRegion();
			Connectivity->SetInputData(Mesh);
			Connectivity->Update();
			Mesh->ShallowCopy(Connectivity->GetOutput());
			Connectivity->Delete();*/
		}

		if (Helper->FillHoles!=0)
		{
			vtkFillHolesFilter *HolesFill= vtkFillHolesFilter::New();
			HolesFill->SetInputData(Mesh);
			HolesFill->SetHoleSize(1e9);
			HolesFill->Update();
			if (HolesFill->GetOutput()->GetNumberOfCells()>0)
			{
				Mesh->ShallowCopy(HolesFill->GetOutput());
			}
			HolesFill->Delete();
		}

		Timer->StopTimer();

		Helper->ProcessingTimes[MyId][0]+=Timer->GetElapsedTime();
		if (Helper->MaximumTimes[MyId][0]<Timer->GetElapsedTime())
			Helper->MaximumTimes[MyId][0]=Timer->GetElapsedTime();

		std::stringstream Name;
		if (Helper->NamesMap[Level].length()>0)
		{
			Name<<Helper->OutputDirectory<<Level<<"-"<<Helper->NamesMap[Level]<<"."<<Helper->GetOutputFormat();
		}
		else
		{
			Name<<Helper->OutputDirectory<<Level<<"."<<Helper->GetOutputFormat();
		}
		int MaxNumberOfVertices=Helper->MaximumNumberOfVertices;
		if ((Mesh->GetNumberOfPoints()>MaxNumberOfVertices)&&(MaxNumberOfVertices>0))
		{
			int WantedNumberOfVertices=MaxNumberOfVertices;
			Timer->StartTimer();
			if (Helper->SimplificationType==1)
			{
				vtkQuadricClustering *Simplification=vtkQuadricClustering::New();
				Simplification->SetInputData(Mesh);
				int NumberOfSubdivisions=(int) pow((double) WantedNumberOfVertices, (double) 1.0/3.0);
				Simplification->SetNumberOfDivisions (NumberOfSubdivisions,NumberOfSubdivisions,NumberOfSubdivisions);
				Simplification->Update();
				Mesh->ShallowCopy(Simplification->GetOutput());
				Simplification->Delete();
			}
			else
			{
				int WantedNumberOfIsotropicVertices=WantedNumberOfVertices;
				if (Helper->Anisotropy!=0)
				{
					WantedNumberOfIsotropicVertices=30*WantedNumberOfVertices;
				}
					
				vtkIsotropicDiscreteRemeshing *Remesh=vtkIsotropicDiscreteRemeshing::New();
				Remesh->GetMetric()->SetGradation(Helper->Gradation);
				vtkSurface *Mesh2=vtkSurface::New();
				Mesh2->CreateFromPolyData(Mesh);
				Timer->StopTimer();
				Helper->ProcessingTimes[MyId][1]+=Timer->GetElapsedTime();
				if (Helper->MaximumTimes[MyId][1]<Timer->GetElapsedTime())
					Helper->MaximumTimes[MyId][1]=Timer->GetElapsedTime();

				Timer->StartTimer();
				Remesh->SetInput(Mesh2);
				if (Helper->ForceManifold!=0)
				{
					Remesh->SetForceManifold(1);
					Remesh->SetSpareFactor(4);
				}
				Remesh->SetNumberOfClusters(WantedNumberOfIsotropicVertices);
				Remesh->SetConsoleOutput(0);
				Remesh->Remesh();

				if (1)
				{
					//Optimize vertices positions with quadrics-based placement
					// Note : this is an adaptation of Siggraph 2000 Paper : Out-of-core simplification of large polygonal models
					vtkIntArray *Clustering=Remesh->GetClustering();
					int Cluster,NumberOfMisclassedItems=0;

					double **ClustersQuadrics =new double*[WantedNumberOfIsotropicVertices];
					for (int i = 0; i < WantedNumberOfIsotropicVertices; i++)
					{
						ClustersQuadrics[i]=new double[9];
						for (int j=0;j<9;j++)
							ClustersQuadrics[i][j]=0;
					}

					vtkIdList *FList=vtkIdList::New();

					for (int i = 0; i < Remesh->GetNumberOfItems (); i++)
					{
						Cluster = Clustering->GetValue (i);
						if ((Cluster >= 0)&& (Cluster < WantedNumberOfIsotropicVertices))
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
					double P[3];
					for (int i = 0; i < WantedNumberOfIsotropicVertices; i++)
					{
						Remesh->GetOutput()->GetPoint (i, P);
						vtkQuadricTools::ComputeRepresentativePoint(ClustersQuadrics[i], P,1);
						Remesh->GetOutput()->SetPointCoordinates (i, P);
						delete[] ClustersQuadrics[i];
					}
					delete [] ClustersQuadrics;
					Mesh->GetPoints()->Modified ();
				}

				if (Helper->Anisotropy!=0)
				{
					vtkAnisotropicDiscreteRemeshing *AnisoRemesh=vtkAnisotropicDiscreteRemeshing::New();
					AnisoRemesh->GetMetric()->SetGradation(Helper->Gradation);
					AnisoRemesh->SetInput(Remesh->GetOutput());
					if (Helper->ForceManifold!=0)
					{
						AnisoRemesh->SetForceManifold(1);
						AnisoRemesh->SetSpareFactor(4);
					}
					AnisoRemesh->SetNumberOfClusters(WantedNumberOfVertices);
					AnisoRemesh->SetConsoleOutput(0);
					AnisoRemesh->Remesh();
					Mesh->ShallowCopy(AnisoRemesh->GetOutput());
					Remesh->Delete();
					AnisoRemesh->Delete();
					Mesh2->Delete();
				}
				else
				{
					Mesh->ShallowCopy(Remesh->GetOutput());
					Remesh->Delete();
					Mesh2->Delete();
				}
			}
			Timer->StopTimer();
			Helper->ProcessingTimes[MyId][2]+=Timer->GetElapsedTime();
			if (Helper->MaximumTimes[MyId][2]<Timer->GetElapsedTime())
				Helper->MaximumTimes[MyId][2]=Timer->GetElapsedTime();
		}

		Timer->StartTimer();

		if (Helper->NumberOfSmoothingSteps!=0)
		{
			vtkWindowedSincPolyDataFilter *Smoother=vtkWindowedSincPolyDataFilter::New();
			Smoother->SetInputData(Mesh);
			Smoother->SetNumberOfIterations(Helper->NumberOfSmoothingSteps);
			Smoother->Update();
			Mesh->ShallowCopy(Smoother->GetOutput());
			Smoother->Delete();
		}

		Writer->SetInputData(Mesh);
		Writer->SetFileName(Name.str().c_str());
		Writer->Write();
		Timer->StopTimer();

		Helper->ProcessingTimes[MyId][4]+=Timer->GetElapsedTime();
		if (Helper->MaximumTimes[MyId][4]<Timer->GetElapsedTime())
			Helper->MaximumTimes[MyId][4]=Timer->GetElapsedTime();

		double GlobalTime=Timer->GetUniversalTime()-StartTime;
		Helper->ProcessingTimes[MyId][5]+=GlobalTime;
		if (Helper->MaximumTimes[MyId][5]<GlobalTime)
			Helper->MaximumTimes[MyId][5]=GlobalTime;		
	}

	Contour->Delete();
	Writer->Delete();
	Image->Delete();
	Timer->Delete();
	return (VTK_THREAD_RETURN_VALUE);
}

template <class voxel_type>
vtkIdList* GetIds(vtkImageData *Image)
{
	vtkIdList *Ids=vtkIdList::New();
	int Dimensions[3];
	Image->GetDimensions(Dimensions);
	int NumVoxels=Dimensions[0]*Dimensions[1]*Dimensions[2];

	voxel_type *Pointer=(voxel_type *) Image->GetScalarPointer();
	voxel_type Label;
	for (;NumVoxels!=0;NumVoxels--)
	{
		 Label=*Pointer;
		 Pointer++;
		if (Ids->IsId(Label)==-1)
			Ids->InsertNextId(Label);
	}
	vtkSortDataArray::Sort(Ids);
	return (Ids);
}

int main( int argc, char *argv[] )
{
	MyThreaderHelperClass Helper;
	vtkMultiThreader *Threader=vtkMultiThreader::New();
	Threader->SetSingleMethod (ThreadedSurfaceExtraction, (void *) &Helper);
	double MaxNumberOfVoxels=100000000;
	int MaximumSubsamplingFactor=4;

	std::map<int, std::string> ColorMap;

	if (argc<2)
	{
		cout<<"Usage : VolumeAnalysis file.mhd [options]"<<endl;
		cout<<"Options : "<<endl;
		cout<<"-n number : sets the number of simplified meshes vertices (default: 32000)"<<endl;
		cout<<"-v number : sets the maximum number of voxels before the image gets resampled (default : 100000000)"<<endl;
		cout<<"-r number : sets the maximum subsampling factor (default : 4)"<<endl;
		cout<<"-j number : sets the number of threads (default : auto)"<<endl;
		cout<<"-s number : sets the simplification method (default: 0 : ACVD)"<<endl;
		cout<<"-sm number : sets the number of smoothing steps (default : 0)"<<endl;
		cout<<"-g number : sets the gradation parameter when using ACVD"<<endl;
		cout<<"-o directory : sets the output directory"<<endl;
		cout<<"-f format : sets the output mesh format (default : vtk)"<<endl;
		cout<<"-m 0/1 : forces manifold output (default : 0)"<<endl;
		cout<<"-a 0/1 : uses anisotropic coarsening (default : 0)"<<endl;
		cout<<"-c 0/1 : keeps (or not) only each biggest component (default : 0)"<<endl;
		cout<<"-x xmlfile : adds colors provided in xmlfile"<<endl;
		exit(1);
	}

	// Parse optionnal arguments
	int ArgumentsIndex=2;
	while (ArgumentsIndex<argc)
	{
		if (strcmp(argv[ArgumentsIndex],"-n")==0)
		{
			cout<<"Setting maximum number of vertices to "<<atoi(argv[ArgumentsIndex+1])<<endl;
			Helper.MaximumNumberOfVertices=(atoi(argv[ArgumentsIndex+1]));
		}
		if (strcmp(argv[ArgumentsIndex],"-j")==0)
		{
			cout<<"Setting number of threads to "<<atoi(argv[ArgumentsIndex+1])<<endl;
			Threader->SetNumberOfThreads(atoi(argv[ArgumentsIndex+1]));
		}
		if (strcmp(argv[ArgumentsIndex],"-s")==0)
		{
			cout<<"Setting simplification type to "<<atoi(argv[ArgumentsIndex+1])<<endl;
			Helper.SimplificationType=atoi(argv[ArgumentsIndex+1]);
		}
		if (strcmp(argv[ArgumentsIndex],"-sm")==0)
		{
			cout<<"Setting number of smoothing steps to "<<atoi(argv[ArgumentsIndex+1])<<endl;
			Helper.NumberOfSmoothingSteps=atoi(argv[ArgumentsIndex+1]);
		}
		if (strcmp(argv[ArgumentsIndex],"-r")==0)
		{
			cout<<"Setting maximum subsampling factor to "<<atoi(argv[ArgumentsIndex+1])<<endl;
			MaximumSubsamplingFactor=atoi(argv[ArgumentsIndex+1]);
		}
		if (strcmp(argv[ArgumentsIndex],"-v")==0)
		{
			cout<<"Setting maximum number of voxels to "<<atoi(argv[ArgumentsIndex+1])<<endl;
			MaxNumberOfVoxels=atoi(argv[ArgumentsIndex+1]);
		}
		if (strcmp(argv[ArgumentsIndex],"-g")==0)
		{
			cout<<"Setting gradation to "<<atof(argv[ArgumentsIndex+1])<<endl;
			Helper.Gradation=atof(argv[ArgumentsIndex+1]);
		}
		if (strcmp(argv[ArgumentsIndex],"-f")==0)
		{
			cout<<"Setting output format to "<<argv[ArgumentsIndex+1]<<endl;
			Helper.SetOutputFormat(argv[ArgumentsIndex+1]);
		}
		if (strcmp(argv[ArgumentsIndex],"-o")==0)
		{
			cout<<"Setting output directory to : "<<argv[ArgumentsIndex+1]<<endl;
			strcpy (Helper.OutputDirectory,argv[ArgumentsIndex+1]);
			strcat (Helper.OutputDirectory,"/");
		}
		if (strcmp(argv[ArgumentsIndex],"-m")==0)
		{
			cout<<"Setting manifold output to "<<atoi(argv[ArgumentsIndex+1])<<endl;
			Helper.ForceManifold=atoi(argv[ArgumentsIndex+1]);
		}
		if (strcmp(argv[ArgumentsIndex],"-fill")==0)
		{
			cout<<"Setting fill holes output to "<<atoi(argv[ArgumentsIndex+1])<<endl;
			Helper.FillHoles=atoi(argv[ArgumentsIndex+1]);
		}
		if (strcmp(argv[ArgumentsIndex],"-a")==0)
		{
			cout<<"Setting anisotropy to "<<atoi(argv[ArgumentsIndex+1])<<endl;
			Helper.Anisotropy=atoi(argv[ArgumentsIndex+1]);
		}
		if (strcmp(argv[ArgumentsIndex],"-c")==0)
		{
			cout<<"Setting keeping biggest component to "<<atoi(argv[ArgumentsIndex+1])<<endl;
			Helper.KeepBiggestComponent=atoi(argv[ArgumentsIndex+1]);
		}
		if (strcmp(argv[ArgumentsIndex],"-x")==0)
		{
			cout<<"Setting colors from : "<<argv[ArgumentsIndex+1]<<endl;
			vtkXMLDataElement *RootElement=vtkXMLUtilities::ReadElementFromFile (argv[ArgumentsIndex+1]);
			for (int j=0;j<RootElement->GetNumberOfNestedElements();j++)
			{
				vtkXMLDataElement *XMLColors=RootElement->GetNestedElement(j);
				if (strcmp(XMLColors->GetName(),"colors")==0)
				{
					for (int i=0;i<XMLColors->GetNumberOfNestedElements();i++)
					{
						vtkXMLDataElement *Element=XMLColors->GetNestedElement(i);
						if (strcmp(Element->GetName(),"color")==0)
						{
							if ((Element->GetAttribute("label")!=0)&&
								(Element->GetAttribute("meshcolor")!=0))
							{
								int Label;
								Element->GetScalarAttribute("label",Label);
								ColorMap[Label]=Element->GetAttribute("meshcolor");
								cout<<"Label "<<Label<<" has color ["<<Element->GetAttribute("meshcolor")<<"]"<<endl;
							}
							if ((Element->GetAttribute("label")!=0)&&
								(Element->GetAttribute("name")!=0))
							{
								int Label;
								Element->GetScalarAttribute("label",Label);
								Helper.NamesMap[Label]=Element->GetAttribute("name");
								cout<<"Label "<<Label<<" has name : "<<Element->GetAttribute("name")<<endl;
							}
						}
					}
					cout<<endl;
				}
			}
		}
		ArgumentsIndex+=2;
	}

	// Load Volume
	cout <<"load : "<<argv[1]<<endl;
	vtkImageData *Image;
	vtkMetaImageReader *Reader=vtkMetaImageReader::New();
	Reader->SetFileName(argv[1]);
	Reader->Update();
	Image=Reader->GetOutput();

	int Dimensions[3];
	Image->GetDimensions(Dimensions);
	cout<<"Image dimensions : "<<Dimensions[0]<<" "<<Dimensions[1]<<" "<<Dimensions[2]<<endl;
	double NumberOfVoxels=Dimensions[0]*Dimensions[1]*Dimensions[2];
	cout<<"Volume = "<<NumberOfVoxels<<" voxels"<<endl;

	vtkIdList *Labels;
	vtkImageCast *Cast;
	switch (Image->GetScalarType())
	{
	case VTK_UNSIGNED_CHAR:
		Labels=GetIds<unsigned char>(Image);
		break;
	case VTK_UNSIGNED_SHORT:
		Labels=GetIds<unsigned short>(Image);
		break;
	case VTK_SHORT:
		Labels=GetIds<short>(Image);
		break;
	case VTK_INT:
		Labels=GetIds<int>(Image);
		break;
	case VTK_SIGNED_CHAR:
		Labels=GetIds<char>(Image);
		break;
	default:
		Cast=vtkImageCast::New();
		Cast->SetInputData(Image);
		Cast->SetOutputScalarTypeToInt ();
		Cast->Update();
		Image->ShallowCopy(Cast->GetOutput());
		Cast->Delete();
		cout<<"Warning : casting image to ints"<<endl;
		Labels=GetIds<int>(Image);
		break;
	}

	cout<<"Found "<<Labels->GetNumberOfIds()<<" labels"<<endl;

	double Factor=pow (MaxNumberOfVoxels/NumberOfVoxels, 1.0/3.0);
	if ((Factor<1)&&(MaximumSubsamplingFactor>1))
	{
		if (Factor<1.0/MaximumSubsamplingFactor)
			Factor=1.0/MaximumSubsamplingFactor;
		vtkImageResample *Resampler=vtkImageResample::New();
		Resampler->SetInputData(Image);
		Resampler->SetInterpolationModeToNearestNeighbor();
		for (int i=0;i<3;i++)
			Resampler->SetAxisMagnificationFactor(i, Factor);
		cout<<"Resampling...";
		Resampler->Update();
		Image=Resampler->GetOutput();
		cout<<"done."<<endl;
		Image->GetDimensions(Dimensions);
		cout<<"New image dimensions : "<<Dimensions[0]<<" "<<Dimensions[1]<<" "<<Dimensions[2]<<endl;
	}

	Helper.Image=Image;
	Helper.Labels=Labels;

	Threader->SingleMethodExecute ();

// Display timings
	cout<<endl<<"Timings :"<<endl;
	for (int TimingType=0;TimingType<NumberOfTimingTypes;TimingType++)
	{
		switch (TimingType)
		{
		case 0:
			cout<<"* Contouring: ";
			break;
		case 1:
			cout<<"* vtkSurface construction: ";
			break;
		case 2:
			cout<<"* Simplification: ";
			break;
		case 3:
			cout<<"* Normal Computation: ";
			break;
		case 4:
			cout<<"* File Output: ";
			break;
		case 5:
		default:
			cout<<"* Whole pipeline: ";
			break;
		}
		double Time=0;
		double MaxTime=0;
		for (int i=0;i<Threader->GetNumberOfThreads();i++)
		{
			double ThreadTime=Helper.ProcessingTimes[i][TimingType];
			Time+=ThreadTime;
			if (MaxTime<Helper.MaximumTimes[i][TimingType])
				MaxTime=Helper.MaximumTimes[i][TimingType];
		}
		cout<<"Average : "<<Time/Labels->GetNumberOfIds()<<" ";
		cout<<"Max : "<<MaxTime<<endl;
	}

	vtkXMLDataElement *Root=vtkXMLDataElement::New ();
	Root->SetName("root");
	Root->SetIntAttribute ("timestamp",(int) vtkTimerLog::GetUniversalTime());

	for (int i=0;i<Labels->GetNumberOfIds();i++)
	{
		std::stringstream Name;
		int Label=Labels->GetId(i);

		if (Helper.NamesMap[Label].length()>0)
		{
			Name<<Label<<"-"<<Helper.NamesMap[Label]<<"."<<Helper.GetOutputFormat();
		}
		else
		{
			Name<<Label<<"."<<Helper.GetOutputFormat();
		}

		vtkXMLDataElement *Element=vtkXMLDataElement::New ();
		Element->SetName("mesh");
		Element->SetIntAttribute("label",Label);
		Element->SetAttribute("mesh",Name.str().c_str());

		std::string Color=ColorMap[Label];
		if (Color.length()==0)
			cout<<"Warning : no color was provided for label "<<Label<<endl;
		else
			Element->SetAttribute("color",Color.c_str());

		Root->AddNestedElement(Element);
		Element->Delete();
	}

	vtkXMLUtilities::WriteElementToFile (Root, "meshes.xml");
	Reader->Delete();
	Root->Delete();
	Threader->Delete();
	Labels->Delete();
    return (0);
}
