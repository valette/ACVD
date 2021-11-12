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

#include <vtkCellData.h>
#include <vtkDiscreteMarchingCubes.h>
#include <vtkFillHolesFilter.h>
#include <vtkIdList.h>
#include <vtkImageCast.h>
#include <vtkImageData.h>
#include <vtkMetaImageWriter.h>
#include <vtkImageExtractComponents.h>
#include <vtkImageResample.h>
#include <vtkImageThreshold.h>
#include <vtkMultiThreader.h>
#include <mutex>
#include <vtkPLYWriter.h>
#include <vtkPointData.h>
#include <vtkPolyDataWriter.h>
#include <vtkQuadricClustering.h>
#include <vtkSortDataArray.h>
#include <vtkSTLWriter.h>
#include <vtkTimerLog.h>
#include <vtkWindowedSincPolyDataFilter.h>
#include <vtkXMLDataElement.h>
#include <vtkXMLUtilities.h>

#include "vtkIsotropicDiscreteRemeshing.h"
#include "vtkAnisotropicDiscreteRemeshing.h"
#include "../VolumeProcessing/vtkRobustImageReader.h"

#define NumberOfTimingTypes 6

class MyThreaderHelperClass
{
public:
    vtkImageData *Image;
    vtkIdList *Labels;
    std::mutex *Lock;
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
		this->Lock=new std::mutex;
		this->MaximumNumberOfVertices=(1<<15) -1;
		this->SimplificationType=0;
		this->NumberOfSmoothingSteps=0;
		strcpy (this->OutputDirectory,"");
		this->Gradation=0;
		this->SetOutputFormat((char *)"vtk");
	}

	~MyThreaderHelperClass()
	{
		delete this->Lock;
	}
};

VTK_THREAD_RETURN_TYPE ThreadedSurfaceExtraction (void *arg)
{
	vtkMultiThreader::ThreadInfo *Info = (vtkMultiThreader::ThreadInfo*) arg;
	MyThreaderHelperClass *Helper = (MyThreaderHelperClass *) Info->UserData;
	int MyId = Info->ThreadID;
	int NumberOfThreads=Info->NumberOfThreads;
	vtkIdList *Labels=Helper->Labels;

	if (MyId >= Labels->GetNumberOfIds()) {
		return (VTK_THREAD_RETURN_VALUE);
	}

	vtkTimerLog *Timer=vtkTimerLog::New();

	for (int i = 0; i < NumberOfTimingTypes; i++) {
		Helper->ProcessingTimes[MyId].push_back(0);
		Helper->MaximumTimes[MyId].push_back(0);
	}

	vtkImageData *Image = vtkImageData::New();
	Image->ShallowCopy(Helper->Image);

	// create the polydataWriter in correct format
	vtkPolyDataWriter *Writer = 0;
	char *Format=Helper->GetOutputFormat();
	char extension[4];
	strcpy (extension,"vtk");
	if (strstr(Format,extension) != NULL) {
		Writer=vtkPolyDataWriter::New();
	} else {
		strcpy (extension, "ply");
		if (strstr(Format,extension) != NULL) {
			Writer = (vtkPolyDataWriter*) vtkPLYWriter::New();
		} else {
			strcpy (extension, "stl");
			if (strstr(Format, extension) != NULL)
				Writer = (vtkPolyDataWriter*) vtkSTLWriter::New();		
		}
	}

	vtkDiscreteMarchingCubes *Contour = vtkDiscreteMarchingCubes::New();
	Contour->ComputeNormalsOff();
	Contour->ComputeGradientsOff ();
	Contour->ComputeScalarsOff();
	Contour->SetInputData(Image);
	for (int i = MyId; i < Labels->GetNumberOfIds(); i += NumberOfThreads) {
		double StartTime = Timer->GetUniversalTime();
		vtkIdType Level = Labels->GetId(i);
		Timer->StartTimer();
		Contour->SetValue(0,Level);
		Contour->Update();
		vtkPolyData *Mesh = Contour->GetOutput();

		if (Helper->KeepBiggestComponent != 0) {
			vtkSurface *ToClean = vtkSurface::New();
			ToClean->CreateFromPolyData(Mesh);
			vtkSurface *Cleaned = ToClean->GetBiggestConnectedComponent();
			Cleaned->EnsureOutwardsNormals();
			Helper->Lock->lock();
			Mesh->ShallowCopy(Cleaned);
			Cleaned->Delete();
			ToClean->Delete();
			Helper->Lock->unlock();
			/*
			vtkPolyDataConnectivityFilter *Connectivity=vtkPolyDataConnectivityFilter::New();
			Connectivity->SetExtractionModeToLargestRegion();
			Connectivity->SetInputData(Mesh);
			Connectivity->Update();
			Mesh->ShallowCopy(Connectivity->GetOutput());
			Connectivity->Delete();*/
		}

		if (Helper->FillHoles != 0) {
			vtkFillHolesFilter *HolesFill = vtkFillHolesFilter::New();
			HolesFill->SetInputData(Mesh);
			HolesFill->SetHoleSize(1e9);
			HolesFill->Update();
			if (HolesFill->GetOutput()->GetNumberOfCells() > 0) {
				Mesh->ShallowCopy(HolesFill->GetOutput());
			}
			HolesFill->Delete();
		}

		Timer->StopTimer();

		Helper->ProcessingTimes[MyId][0] += Timer->GetElapsedTime();
		if (Helper->MaximumTimes[MyId][0] < Timer->GetElapsedTime())
			Helper->MaximumTimes[MyId][0] = Timer->GetElapsedTime();

		std::stringstream Name;
		Name << Helper->OutputDirectory << Level;
		if (Helper->NamesMap[Level].length()>0) {
			Name << "-" << Helper->NamesMap[Level] << "." << Helper->GetOutputFormat();
		} else {
			Name << "." << Helper->GetOutputFormat();
		}
		int MaxNumberOfVertices = Helper->MaximumNumberOfVertices;
		if ((Mesh->GetNumberOfPoints() > MaxNumberOfVertices) && (MaxNumberOfVertices > 0))
		{
			int WantedNumberOfVertices = MaxNumberOfVertices;
			Timer->StartTimer();
			if (Helper->SimplificationType == 1) {
				vtkQuadricClustering *Simplification = vtkQuadricClustering::New();
				Simplification->SetInputData(Mesh);
				int NumberOfSubdivisions=(int) pow((double) WantedNumberOfVertices, (double) 1.0/3.0);
				Simplification->SetNumberOfDivisions (NumberOfSubdivisions,NumberOfSubdivisions,NumberOfSubdivisions);
				Simplification->Update();
				Mesh->ShallowCopy(Simplification->GetOutput());
				Simplification->Delete();
			} else {
				int WantedNumberOfIsotropicVertices=WantedNumberOfVertices;
				if (Helper->Anisotropy!=0) {
					WantedNumberOfIsotropicVertices=30*WantedNumberOfVertices;
				}
					
				vtkIsotropicDiscreteRemeshing *Remesh=vtkIsotropicDiscreteRemeshing::New();
				Remesh->GetMetric()->SetGradation(Helper->Gradation);
				vtkSurface *Mesh2 = vtkSurface::New();
				Mesh2->CreateFromPolyData(Mesh);
				Timer->StopTimer();
				Helper->ProcessingTimes[MyId][1] += Timer->GetElapsedTime();
				if (Helper->MaximumTimes[MyId][1] < Timer->GetElapsedTime())
					Helper->MaximumTimes[MyId][1] = Timer->GetElapsedTime();

				Timer->StartTimer();
				Remesh->SetInput(Mesh2);
				if (Helper->ForceManifold!=0) {
					Remesh->SetForceManifold(1);
				}
				Remesh->SetNumberOfClusters(WantedNumberOfIsotropicVertices);
				Remesh->SetConsoleOutput(0);
				Remesh->SetSubsamplingThreshold(1);
				Remesh->Remesh();

				if (1) {
					//Optimize vertices positions with quadrics-based placement
					// Note : this is an adaptation of Siggraph 2000 Paper : Out-of-core simplification of large polygonal models
					vtkIntArray *Clustering=Remesh->GetClustering();
					int Cluster,NumberOfMisclassedItems=0;

					double **ClustersQuadrics =new double*[WantedNumberOfIsotropicVertices];
					for (int i = 0; i < WantedNumberOfIsotropicVertices; i++) {
						ClustersQuadrics[i] = new double[9];
						for (int j = 0; j < 9; j++) {
							ClustersQuadrics[i][j] = 0;
						}
					}

					vtkIdList *FList = vtkIdList::New();

					for (int i = 0; i < Remesh->GetNumberOfItems (); i++) {
						Cluster = Clustering->GetValue (i);
						if ((Cluster >= 0)&& (Cluster < WantedNumberOfIsotropicVertices)) {
							if (Remesh->GetClusteringType() == 0) {
								vtkQuadricTools::AddTriangleQuadric(ClustersQuadrics[Cluster],Remesh->GetInput(),i,false);
							} else {
								Remesh->GetInput()->GetVertexNeighbourFaces(i,FList);
								for (int j=0;j<FList->GetNumberOfIds();j++) {
									vtkQuadricTools::AddTriangleQuadric(ClustersQuadrics[Cluster]
											,Remesh->GetInput(),FList->GetId(j),false);
								}
							}
						} else {
							NumberOfMisclassedItems++;
						}
					}
					FList->Delete();
					double P[3];
					for (int i = 0; i < WantedNumberOfIsotropicVertices; i++) {
						Remesh->GetOutput()->GetPoint (i, P);
						vtkQuadricTools::ComputeRepresentativePoint(ClustersQuadrics[i], P,1);
						Remesh->GetOutput()->SetPointCoordinates (i, P);
						delete[] ClustersQuadrics[i];
					}
					delete [] ClustersQuadrics;
					Mesh->GetPoints()->Modified ();
				}

				if (Helper->Anisotropy != 0) {
					vtkAnisotropicDiscreteRemeshing *AnisoRemesh=vtkAnisotropicDiscreteRemeshing::New();
					AnisoRemesh->GetMetric()->SetGradation(Helper->Gradation);
					AnisoRemesh->SetInput(Remesh->GetOutput());
					if (Helper->ForceManifold!=0) {
						AnisoRemesh->SetForceManifold(1);
					}
					AnisoRemesh->SetNumberOfClusters(WantedNumberOfVertices);
					AnisoRemesh->SetConsoleOutput(0);
					AnisoRemesh->Remesh();
					Helper->Lock->lock();
					Mesh->ShallowCopy(AnisoRemesh->GetOutput());
					Remesh->Delete();
					AnisoRemesh->Delete();
					Mesh2->Delete();
					Helper->Lock->unlock();
				} else {
					Helper->Lock->lock();
					Mesh->ShallowCopy(Remesh->GetOutput());
					Remesh->Delete();
					Mesh2->Delete();
					Helper->Lock->unlock();
				}
			}
			Timer->StopTimer();
			Helper->ProcessingTimes[MyId][2]+=Timer->GetElapsedTime();
			if (Helper->MaximumTimes[MyId][2]<Timer->GetElapsedTime())
				Helper->MaximumTimes[MyId][2]=Timer->GetElapsedTime();
		}

		Timer->StartTimer();

		if (Helper->NumberOfSmoothingSteps != 0) {
			vtkWindowedSincPolyDataFilter *Smoother=vtkWindowedSincPolyDataFilter::New();
			Smoother->SetInputData(Mesh);
			Smoother->SetNumberOfIterations(Helper->NumberOfSmoothingSteps);
			Smoother->Update();
			Mesh->ShallowCopy(Smoother->GetOutput());
			Helper->Lock->lock();
			Smoother->Delete();
			Helper->Lock->unlock();
		}

		Writer->SetInputData(Mesh);
		Writer->SetFileName(Name.str().c_str());
		Writer->Write();
		Timer->StopTimer();
		
		Helper->Lock->lock();
		cout<< "Label " << Level << " done"<<endl;
		Helper->Lock->unlock();

		Helper->ProcessingTimes[MyId][4]+=Timer->GetElapsedTime();
		if (Helper->MaximumTimes[MyId][4]<Timer->GetElapsedTime())
			Helper->MaximumTimes[MyId][4]=Timer->GetElapsedTime();

		double GlobalTime=Timer->GetUniversalTime()-StartTime;
		Helper->ProcessingTimes[MyId][5]+=GlobalTime;
		if (Helper->MaximumTimes[MyId][5]<GlobalTime)
			Helper->MaximumTimes[MyId][5]=GlobalTime;		
	}

	Helper->Lock->lock();
	Contour->Delete();
	Writer->Delete();
	Image->Delete();
	Timer->Delete();
	cout<< "Thread " << MyId << " done"<<endl;
	Helper->Lock->unlock();
	return (VTK_THREAD_RETURN_VALUE);
}

template <class voxel_type>
vtkIdList* GetIds(vtkImageData *Image)
{
	vtkIdList *Ids=vtkIdList::New();
	int Dimensions[3];
	Image->GetDimensions(Dimensions);
	int NumVoxels = Dimensions[0] * Dimensions[1] * Dimensions[2];

	voxel_type *Pointer=(voxel_type *) Image->GetScalarPointer();
	voxel_type Label;
	for (;NumVoxels!=0;NumVoxels--) {
		 Label = *Pointer;
		 Pointer++;
		if (Ids->IsId(Label) == -1) {
			Ids->InsertNextId(Label);
		}
	}
	vtkSortDataArray::Sort(Ids);
	return (Ids);
}

int main( int argc, char *argv[] )
{
	MyThreaderHelperClass Helper;
	vtkMultiThreader *Threader=vtkMultiThreader::New();
	Threader->SetSingleMethod (ThreadedSurfaceExtraction, (void *) &Helper);
	double MaxNumberOfVoxels = VTK_INT_MAX;
	int MaximumSubsamplingFactor=4;
	double threshold = VTK_DOUBLE_MIN;
	int binaryMasks = 0;

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
		cout<<"-t threshold_value : apply thresholding to the volume (default : not used)"<<endl;
		cout<<"-x xmlfile : adds colors provided in xmlfile"<<endl;
		exit(1);
	}

	// Parse optionnal arguments
	int ArgumentsIndex=2;
	while (ArgumentsIndex < argc) {
		char *key = argv[ArgumentsIndex];
		char *value = argv[ArgumentsIndex + 1];

		if (strcmp(key, "-masks") == 0) {
			cout << "Output binary masks : " << atoi(value) << endl;
			binaryMasks = atoi(value);
		}
		if (strcmp(key, "-n") == 0) {
			cout << "Setting maximum number of vertices to " << atoi(value) << endl;
			Helper.MaximumNumberOfVertices = (atoi(value));
		}
		if (strcmp(key, "-j") == 0) {
			cout << "Setting number of threads to " << atoi(value) << endl;
			Threader->SetNumberOfThreads(atoi(value));
		}
		if (strcmp(key, "-s") == 0) {
			cout << "Setting simplification type to " << atoi(value) << endl;
			Helper.SimplificationType=atoi(value);
		}
		if (strcmp(key, "-sm") == 0) {
			cout << "Setting number of smoothing steps to " << atoi(value) << endl;
			Helper.NumberOfSmoothingSteps=atoi(value);
		}
		if (strcmp(key, "-r") == 0) {
			cout << "Setting maximum subsampling factor to " << atoi(value) << endl;
			MaximumSubsamplingFactor=atoi(value);
		}
		if (strcmp(key, "-v") == 0) {
			cout << "Setting maximum number of voxels to " << atoi(value) << endl;
			MaxNumberOfVoxels=atoi(value);
		}
		if (strcmp(key, "-g") == 0) {
			cout << "Setting gradation to " << atof(value) << endl;
			Helper.Gradation=atof(value);
		}
		if (strcmp(key, "-f") == 0) {
			cout << "Setting output format to " << value << endl;
			Helper.SetOutputFormat(value);
		}
		if (strcmp(key, "-o") == 0) {
			cout << "Setting output directory to : " << value << endl;
			strcpy (Helper.OutputDirectory,value);
			strcat (Helper.OutputDirectory,"/");
		}
		if (strcmp(key, "-m") == 0) {
			cout<<"Setting manifold output to " << atoi(value) << endl;
			Helper.ForceManifold=atoi(value);
		}
		if (strcmp(key, "-fill") == 0) {
			cout << "Setting fill holes output to " << atoi(value) << endl;
			Helper.FillHoles=atoi(value);
		}
		if (strcmp(key, "-a") == 0) {
			cout << "Setting anisotropy to " << atoi(value) << endl;
			Helper.Anisotropy=atoi(value);
		}
		if (strcmp(key, "-c") == 0) {
			cout << "Setting keeping biggest component to " << atoi(value) << endl;
			Helper.KeepBiggestComponent=atoi(value);
		}
		if (strcmp(key, "-t") == 0) {
			cout << "Threshold value : " << atof(value) << endl;
			threshold = atoi(value);
		}
		if (strcmp(key, "-x") == 0) {
			cout << "Setting colors from : " << value << endl;
			vtkXMLDataElement *RootElement=vtkXMLUtilities::ReadElementFromFile (value);
			for (int j=0;j<RootElement->GetNumberOfNestedElements();j++) {
				vtkXMLDataElement *XMLColors=RootElement->GetNestedElement(j);
				if (strcmp(XMLColors->GetName(),"colors")==0) {
					for (int i=0;i<XMLColors->GetNumberOfNestedElements();i++) {
						vtkXMLDataElement *Element=XMLColors->GetNestedElement(i);
						if (strcmp(Element->GetName(),"color")==0) {
							if ((Element->GetAttribute("label")!=0) &&
								(Element->GetAttribute("meshcolor")!=0)) {
								int Label;
								Element->GetScalarAttribute("label",Label);
								ColorMap[Label]=Element->GetAttribute("meshcolor");
								cout<<"Label "<<Label<<" has color ["<<Element->GetAttribute("meshcolor")<<"]"<<endl;
							}
							if ((Element->GetAttribute("label")!=0) &&
								(Element->GetAttribute("name")!=0)) {
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
	vtkRobustImageReader *Reader = vtkRobustImageReader::New();
	Reader->SetFileName(argv[1]);
	Reader->Update();
	Image = Reader->GetOutput();
	int numComp = Image->GetNumberOfScalarComponents();

	if ( numComp > 1) {
		cout << "Warning : image has " << numComp << " components. Will keep only the first" << endl;
		vtkImageExtractComponents *components = vtkImageExtractComponents::New();
		components->SetInputData(Image);
		components->SetComponents(0);
		components->Update();
		Image = components->GetOutput();
	}

	if (threshold > VTK_DOUBLE_MIN) {
		vtkImageThreshold *t = vtkImageThreshold::New();
		t->SetOutputScalarTypeToUnsignedChar();
		t->SetInValue(1);
		t->SetOutValue(0);
		t->SetInputData(Image);
		t->ThresholdByUpper(threshold);
		t->Update();
		Image = t->GetOutput();
	}

	int Dimensions[3];
	Image->GetDimensions(Dimensions);
	cout<<"Image dimensions : "<<Dimensions[0]<<" "<<Dimensions[1]<<" "<<Dimensions[2]<<endl;
	double NumberOfVoxels=Dimensions[0]*Dimensions[1]*Dimensions[2];
	cout<<"Volume = "<<NumberOfVoxels<<" voxels"<<endl;

	vtkIdList *Labels;
	vtkImageCast *Cast;
	switch (Image->GetScalarType()) {
	case VTK_UNSIGNED_CHAR:
		Labels = GetIds<unsigned char>(Image);
		break;
	case VTK_UNSIGNED_SHORT:
		Labels = GetIds<unsigned short>(Image);
		break;
	case VTK_SHORT:
		Labels = GetIds<short>(Image);
		break;
	case VTK_INT:
		Labels = GetIds<int>(Image);
		break;
	case VTK_SIGNED_CHAR:
		Labels = GetIds<char>(Image);
		break;
	default:
		Cast = vtkImageCast::New();
		Cast->SetInputData(Image);
		Cast->SetOutputScalarTypeToInt ();
		Cast->Update();
		Image->ShallowCopy(Cast->GetOutput());
		Cast->Delete();
		cout << "Warning : casting image to ints" << endl;
		Labels = GetIds<int>(Image);
		break;
	}

	cout << "Found " << Labels->GetNumberOfIds() << " labels" << endl;

	if ( binaryMasks ) {

		vtkImageThreshold *threshold = vtkImageThreshold::New();
		threshold->SetInValue(1);
		threshold->SetOutValue(0);
		threshold->SetOutputScalarTypeToUnsignedChar();
		threshold->SetInputData( Image );
		vtkMetaImageWriter *writer = vtkMetaImageWriter::New();

		for (int i = 0; i < Labels->GetNumberOfIds(); i++) {

			std::stringstream name;
			int label = Labels->GetId( i );
			name << "label" << label << ".mhd";
			threshold->ThresholdBetween( label, label );
			threshold->Update();
			writer->SetInputData( threshold->GetOutput() );
			writer->SetFileName( name.str().c_str() );
			writer->Write();

		}

	}


	double Factor = pow (MaxNumberOfVoxels / NumberOfVoxels, 1.0 / 3.0);
	if ((Factor < 1) && (MaximumSubsamplingFactor > 1))
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
	Helper.ProcessingTimes.resize(Threader->GetNumberOfThreads());
	Helper.MaximumTimes.resize(Threader->GetNumberOfThreads());
	Threader->SingleMethodExecute ();
	cout << "All meshes extracted" << endl;
// Display timings
	cout<<endl<<"Timings :"<<endl;
	for (int TimingType=0;TimingType<NumberOfTimingTypes;TimingType++) {
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
		for (int i=0;i<Threader->GetNumberOfThreads();i++) {
			if (Helper.ProcessingTimes[i].size() <= i) continue;
			double ThreadTime = Helper.ProcessingTimes[i][TimingType];
			Time += ThreadTime;
			if (MaxTime < Helper.MaximumTimes[i][TimingType])
				MaxTime = Helper.MaximumTimes[i][TimingType];
		}
		cout<<"Average : " << Time/Labels->GetNumberOfIds() << " ";
		cout<<"Max : " << MaxTime << endl;
	}

	vtkXMLDataElement *Root=vtkXMLDataElement::New ();
	Root->SetName("root");
	Root->SetIntAttribute ("timestamp",(int) vtkTimerLog::GetUniversalTime());

	for (int i = 0; i < Labels->GetNumberOfIds(); i++) {
		std::stringstream Name;
		int Label = Labels->GetId(i);

		if (Helper.NamesMap[Label].length() > 0) {
			Name << Label<< "-" << Helper.NamesMap[Label]<<"."<<Helper.GetOutputFormat();
		} else {
			Name << Label << "."<< Helper.GetOutputFormat();
		}

		vtkXMLDataElement *Element=vtkXMLDataElement::New ();
		Element->SetName("mesh");
		Element->SetIntAttribute("label",Label);
		Element->SetAttribute("mesh",Name.str().c_str());

		std::string Color=ColorMap[Label];
		if (Color.length()==0) {
			cout<<"Warning : no color was provided for label "<<Label<<endl;
		} else {
			Element->SetAttribute("color",Color.c_str());
		}

		Root->AddNestedElement(Element);
		Element->Delete();
		cout<< "XML element for label " << Label << " done"<<endl;
	}

	vtkXMLUtilities::WriteElementToFile (Root, "meshes.xml");
//	Reader->Delete();
	cout << "reader deleted" << endl;
//	Root->Delete();
	cout << "root deleted" << endl;
//	Threader->Delete();
	cout << "threader deleted" << endl;
//	Labels->Delete();
	cout << "labels deleted" << endl;
    return (0);
}
