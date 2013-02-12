/*=========================================================================

Program:   VolumeSlice : generates images slices from volume data
Module:    vtkSurface
Language:  C++
Date:      2010/03
Auteur:   Sebastien Valette

=========================================================================*/
// .NAME VolumeSlice
// .SECTION Description

#include <sstream>
#include <cstring>

#include <cstdlib>
#include <string>
#include <fstream>
#include <iostream>

#include <vtkPNGWriter.h>
#include <vtkJPEGWriter.h>
#include <vtkImageData.h>
#include <vtkMetaImageReader.h>
#include <vtkImageReslice.h>
#include <vtkImageShiftScale.h>
#include <vtkLookupTable.h>
#include <vtkImageMapToColors.h>
#include <vtkTIFFWriter.h>
#include <vtkXMLDataElement.h>
#include <vtkXMLUtilities.h>
#include <vtkMultiThreader.h>
#include <vtkTimerLog.h>

using std::cout ;
using std::cerr ;
using std::string ;
using std::ostream ;
using std::fstream ;
using std::ifstream ;
using std::ofstream ;
using std::istringstream ;
using std::ostringstream ;

class MyThreaderHelperClass
{
public:
    vtkImageData *Image;

	char prefix[1000];
	char directory[1000];
	int offset;
	double Scale;
	double Shift;
	bool userdefinedscale;
	int Format;
	double Range[2];
	int NumberOfSlices;
	bool pngUsingMHDColors;
	int xyzOrientation;
	double	nativeTransformMtrx[3][3];
	char anatomicalOrientation[1000];

	MyThreaderHelperClass()
	{
		strcpy (prefix,"slice");
		strcpy (directory, "");
		offset=10000;
		userdefinedscale=false;
		Shift=0;
		Scale=1;
		Format=1;
		pngUsingMHDColors=false;
		xyzOrientation=0;
	}

	~MyThreaderHelperClass()
	{
	}
};

VTK_THREAD_RETURN_TYPE ThreadedImageSlice (void *arg)
{
	vtkMultiThreader::ThreadInfo *Info = (vtkMultiThreader::ThreadInfo*) arg;
	MyThreaderHelperClass *Helper = (MyThreaderHelperClass *) Info->UserData;
	int MyId = Info->ThreadID;
	int NumberOfThreads=Info->NumberOfThreads;
	
	vtkImageData *Image=vtkImageData::New();
	Image->ShallowCopy(Helper->Image);
	vtkImageReslice *Slicer=vtkImageReslice::New();
	Slicer->SetNumberOfThreads(1);
	vtkImageShiftScale *Cast=vtkImageShiftScale::New();
	Cast->SetNumberOfThreads(1);
	vtkImageWriter *Writer;
	std::stringstream Suffix;

	int Dimensions[3];
	int orientationDims[3];
	int NumberOfSlices;
	int Point[3];
	int pointIndex;
	double params[3][3];
	std::string axesName;
	switch(Helper->xyzOrientation)
	{
	//// Slice on x (x decreasing): ZY X
		case 1 :
			Image->GetDimensions(orientationDims);
			NumberOfSlices = orientationDims[0];
			Dimensions[0] = orientationDims[2];
			Dimensions[1] = orientationDims[1];
			Dimensions[2] = 1;
			params[0][0]=0;		params[0][1]=0;		params[0][2]=1;
			params[1][0]=0;		params[1][1]=-1;	params[1][2]=0;
			params[2][0]=1;		params[2][1]=0;		params[2][2]=0;
			Point[2] = 0;
			Point[1] = 0;
			pointIndex = 0;
			axesName = "ZY";
			break;
	//// Slice on y (y growing): XZ Y
		case 2 :
			Image->GetDimensions(orientationDims);
			NumberOfSlices = orientationDims[1];
			Dimensions[0] = orientationDims[0];
			Dimensions[1] = orientationDims[2];
			Dimensions[2] = 1;
			params[0][0]=1;		params[0][1]=0;		params[0][2]=0;
			params[1][0]=0;		params[1][1]=0;		params[1][2]=-1;
			params[2][0]=0;		params[2][1]=-1;	params[2][2]=0;
			Point[0] = 0;
			Point[2] = 0;
			pointIndex = 1;
			axesName = "XZ";
			break;
	//// Slice on z (z growing): XY Z
		default :
			Image->GetDimensions(Dimensions);
			NumberOfSlices = Dimensions[2];
			Dimensions[2] = 1;
			params[0][0]=1;		params[0][1]=0;		params[0][2]=0;
			params[1][0]=0;		params[1][1]=-1;	params[1][2]=0;
			params[2][0]=0;		params[2][1]=0;		params[2][2]=1;
			Point[0] = 0;
			Point[1] = 0;
			pointIndex = 2;
			axesName = "XY";
	};
	cout<<" "<<endl;
	cout<<">>>>>>> params :"<<endl;
	for (int i=0; i<3; ++i)
	{
		for(int j=0; j<3; ++j)
			cout<<params[i][j]<<",";
		cout<<endl;
	}

	Slicer->SetResliceAxesDirectionCosines(	params[0][0],	params[0][1],	params[0][2],
											params[1][0],	params[1][1],	params[1][2],
											params[2][0],	params[2][1],	params[2][2]);
	vtkImageData *PNGSlice=0;


	int NumComponents;
	switch (Image->GetScalarType())
	{
	case 2:
	case 15:
	case 3:
		NumComponents=1;
		break;
	case 4:
	case 5:
		NumComponents=2;
		break;
	default:
		NumComponents=4;
		break;
	}

	if (Helper->Format==0)
	{
		Cast->SetOutputScalarTypeToFloat ();
		Writer=vtkPNGWriter::New();
		Suffix<<".png";
		PNGSlice=vtkImageData::New();
		PNGSlice->SetNumberOfScalarComponents(NumComponents);
		PNGSlice->SetScalarTypeToUnsignedChar();
		PNGSlice->SetDimensions(Dimensions);
		PNGSlice->AllocateScalars();		
	}
	else
	{
		Cast->SetOutputScalarTypeToUnsignedChar ();
		Writer=vtkJPEGWriter::New();
		Suffix<<".jpg";
		if (!Helper->userdefinedscale)
		{
			Cast->SetShift(-Helper->Range[0]);
			Cast->SetScale(255/(Helper->Range[1]-Helper->Range[0]));
		}
		else
		{
			Cast->SetShift(Helper->Shift);
			Cast->SetScale(Helper->Scale);
		}
	}
	
	Slicer->SetInput(Image);
	Slicer->SetOutputDimensionality(2);

	for (int i=MyId;i<NumberOfSlices;i+=NumberOfThreads)
	{
		Point[pointIndex]=i;
		double Coords[3];
		Image->GetPoint (Image->ComputePointId(Point), Coords);
		Slicer->SetResliceAxesOrigin(Coords);
		Slicer->Update();
		if (Image->GetNumberOfScalarComponents()>1){
			Writer->SetInput(Slicer->GetOutput());
		}
		else
		{
			if (Helper->Format==0)
			{
				vtkImageData *Output;
				unsigned char* SlicePointer;
				int j;
				int numPixels;

				switch (Image->GetScalarType())
				{
				case 2:
				case 15:
				case 3:
				case 4:
				case 5:
					Output=Slicer->GetOutput();
					break;
				default:
					Cast->SetInput(Slicer->GetOutput());
					Cast->Update();
					Output=Cast->GetOutput();
					break;
				}

				unsigned char* PNGPointer=(unsigned char *) PNGSlice->GetScalarPointer();
				SlicePointer=(unsigned char *) Output->GetScalarPointer();
				numPixels=Dimensions[0]*Dimensions[1]*NumComponents;
				for (j=0;j!=numPixels;j++)
				{
					PNGPointer[j]=SlicePointer[j];
				}

				Writer->SetInput(PNGSlice);
			}
			else
			{
				Cast->SetInput(Slicer->GetOutput());
				Cast->Update();
				Writer->SetInput(Cast->GetOutput());
			}
		}
		std::stringstream Name;
		Name<<Helper->directory<<Helper->prefix<<axesName<<i+Helper->offset<<Suffix.str();
		Writer->SetFileName(Name.str().c_str());
		Writer->Write();
	}

	if (Helper->Format==0)
	{
		PNGSlice->Delete();
	}

	return (VTK_THREAD_RETURN_VALUE);
}

int main( int argc, char *argv[] )
{
	MyThreaderHelperClass Helper;

	if (argc<2)
	{
		cout<<"Usage : VolumeSlice file.mhd [options]"<<endl;
		cout<<"Options : "<<endl;
		cout<<"-f value  : set format (0: png ; 1 : jpg) default : 1"<<endl;
		cout<<"-sc value : set scale "<<endl;
		cout<<"-sh value : set shift "<<endl;
		cout<<"-o value : set output directory"<<endl;
		cout<<"-mhdcolors : create png image without using scale nor shifting"<<endl;
		cout<<"-orientation : set the axe used for the orientation of the slice (0 : z; 1 : x; 2 : y) default : 0"<<endl;
		exit(1);
	}
	
	// Load Volume
	cout <<"load : "<<argv[1]<<endl;
	
	vtkMetaImageReader *Reader=vtkMetaImageReader::New();
	Reader->SetFileName(argv[1]);
	Reader->Update();
	vtkImageData *Image=Reader->GetOutput();
	Image->GetScalarRange(Helper.Range);
	cout<<"Range : "<<Helper.Range[0]<<"   "<<Helper.Range[1]<<endl;
	if (Image->GetNumberOfScalarComponents()>1)
		Helper.pngUsingMHDColors = true;
	
	Helper.Image=Image;
	
	bool matrixFound = false;
	bool anOrientation = false;
	string matrixPrefix;
	// Keep transformation matrix and anatomical orientation of input meta image, it needs to be copied to output images.
	// Since VTK currently offers no methods to read and write this information, we're forced to use IO streams
	string imageReaderType = Reader->GetClassName();
	for(int i=0; i<3; ++i)
		for(int j=0; j<3; ++j)
			Helper.nativeTransformMtrx[j][i] = 0;
	Helper.nativeTransformMtrx[0][0] = 1;
	Helper.nativeTransformMtrx[1][1] = 1;
	Helper.nativeTransformMtrx[2][2] = 1;
	if(imageReaderType=="vtkMetaImageReader")
	{
		ifstream is(argv[1]);
		string line, prefix;
		while(is)
		{
			getline(is, line);
			istringstream iss(line);
			iss >> prefix;
			if((prefix=="TransformMatrix")||(prefix == "Orientation")||(prefix == "Rotation"))
			{
				matrixPrefix = prefix;
				iss >> prefix; // eat the "=" sign
				for(int i=0; i<3; ++i)
					for(int j=0; j<3; ++j)
						iss >> Helper.nativeTransformMtrx[j][i];
				matrixFound = true;
			}
			if(prefix=="AnatomicalOrientation")
			{
				iss >> prefix; // eat the "=" sign
				iss >> Helper.anatomicalOrientation;
				anOrientation = true;
			}
		}
		is.close();
	}

	// Parse optionnal arguments
	int ArgumentsIndex=2;
	while (ArgumentsIndex<argc)
	{
		if (strcmp(argv[ArgumentsIndex],"-sc")==0)
		{
			double Scale=atof(argv[ArgumentsIndex+1]);
			cout<<"Scale="<<Scale<<endl;
			Helper.Scale=Scale;
			Helper.userdefinedscale=true;
		}

		if (strcmp(argv[ArgumentsIndex],"-sh")==0)
		{
			double Shift=atof(argv[ArgumentsIndex+1]);
			cout<<"Shift="<<Shift<<endl;
			Helper.Shift=Shift;
			Helper.userdefinedscale=true;
		}

		if (strcmp(argv[ArgumentsIndex],"-f")==0)
		{

			Helper.Format=atoi(argv[ArgumentsIndex+1]);
			cout<<"Format: "<<Helper.Format<<endl;
		}

		if (strcmp(argv[ArgumentsIndex],"-o")==0)
		{
			strcpy (Helper.directory, argv[ArgumentsIndex+1]);
			cout<<"Output directory : "<< Helper.directory <<endl;
		}

		if (strcmp(argv[ArgumentsIndex],"-mhdcolors")==0)
		{
			Helper.pngUsingMHDColors = true;
		}
		
		if (strcmp(argv[ArgumentsIndex],"-orientation")==0)
		{
			int tempVal = atoi(argv[ArgumentsIndex+1]);
			if((0<=tempVal)&&(tempVal<3))
				Helper.xyzOrientation = tempVal;
		}
		
		ArgumentsIndex+=2;
	}

	int Dimensions[3];
	Image->GetDimensions(Dimensions);
	//// Slice on z (z growing):
	//~ Helper.NumberOfSlices=Dimensions[2];
	//// Slice on x (x decreasing):
	//~ Helper.NumberOfSlices=Dimensions[0];
	//// Slice on y (y growing):
	//~ Helper.NumberOfSlices=Dimensions[1];
	
	cout<<"Dimensions : "<<Dimensions[0]<<" "<<Dimensions[1]<<" "<<Dimensions[2]<<endl;

	vtkMultiThreader *Threader=vtkMultiThreader::New();
	Threader->SetSingleMethod (ThreadedImageSlice, (void *) &Helper);
	cout<<"Using "<<Threader->GetNumberOfThreads()<<" threads"<<endl;
	Threader->SingleMethodExecute ();
	Threader->Delete();

	//output xml file containint volume details
	vtkXMLDataElement *Root=vtkXMLDataElement::New ();
	Root->SetName("volume");

	vtkXMLDataElement *OriginalName=vtkXMLDataElement::New ();
	OriginalName->SetName("name");
	OriginalName->SetCharacterData (argv[1], strlen(argv[1])+1);
	Root->AddNestedElement(OriginalName);

	vtkXMLDataElement *XMLDimensions=vtkXMLDataElement::New ();
	XMLDimensions->SetName("dimensions");
	XMLDimensions->SetIntAttribute ("x", Dimensions[0]);
	XMLDimensions->SetIntAttribute ("y", Dimensions[1]);
	XMLDimensions->SetIntAttribute ("z", Dimensions[2]);
	Root->AddNestedElement(XMLDimensions);

	vtkXMLDataElement *XMLRange=vtkXMLDataElement::New ();
	XMLRange->SetName("scalars");
	XMLRange->SetIntAttribute("numberOfScalarComponents", Image->GetNumberOfScalarComponents());
	XMLRange->SetFloatAttribute ("min", Helper.Range[0]);
	XMLRange->SetFloatAttribute ("max", Helper.Range[1]);
	XMLRange->SetIntAttribute ("type", 	Image->GetScalarType());
	XMLRange->SetIntAttribute ("size", 	Image->GetScalarSize());
	char typeString[1000];
	strcpy(typeString,Image->GetScalarTypeAsString());
	XMLRange->SetCharacterData (typeString, strlen(typeString)+1);
	Root->AddNestedElement(XMLRange);

	double Origin[3];
	Image->GetOrigin(Origin);
	vtkXMLDataElement *XMLOrigin=vtkXMLDataElement::New ();
	XMLOrigin->SetName("origin");
	XMLOrigin->SetDoubleAttribute ("x", Origin[0]);
	XMLOrigin->SetDoubleAttribute ("y", Origin[1]);
	XMLOrigin->SetDoubleAttribute ("z", Origin[2]);
	Root->AddNestedElement(XMLOrigin);

	double Spacing[3];
	Image->GetSpacing(Spacing);
	vtkXMLDataElement *XMLSpacing=vtkXMLDataElement::New ();
	XMLSpacing->SetName("spacing");
	XMLSpacing->SetDoubleAttribute ("x", Spacing[0]);
	XMLSpacing->SetDoubleAttribute ("y", Spacing[1]);
	XMLSpacing->SetDoubleAttribute ("z", Spacing[2]);
	Root->AddNestedElement(XMLSpacing);

	int Extent[6];
	Image->GetExtent(Extent);
	vtkXMLDataElement *XMLExtent=vtkXMLDataElement::New ();
	XMLExtent->SetName("extent");
	XMLExtent->SetIntAttribute ("x1", Extent[0]);
	XMLExtent->SetIntAttribute ("x2", Extent[1]);
	XMLExtent->SetIntAttribute ("y1", Extent[2]);
	XMLExtent->SetIntAttribute ("y2", Extent[3]);
	XMLExtent->SetIntAttribute ("z1", Extent[4]);
	XMLExtent->SetIntAttribute ("z2", Extent[5]);
	Root->AddNestedElement(XMLExtent);

	if(matrixFound)
	{
		vtkXMLDataElement *XMLMatrix = vtkXMLDataElement::New ();
		XMLMatrix->SetName("matrix");
		XMLMatrix->SetDoubleAttribute ("x1", Helper.nativeTransformMtrx[0][0]);
		XMLMatrix->SetDoubleAttribute ("x2", Helper.nativeTransformMtrx[1][0]);
		XMLMatrix->SetDoubleAttribute ("x3", Helper.nativeTransformMtrx[2][0]);
		XMLMatrix->SetDoubleAttribute ("y1", Helper.nativeTransformMtrx[0][1]);
		XMLMatrix->SetDoubleAttribute ("y2", Helper.nativeTransformMtrx[1][1]);
		XMLMatrix->SetDoubleAttribute ("y3", Helper.nativeTransformMtrx[2][1]);
		XMLMatrix->SetDoubleAttribute ("z1", Helper.nativeTransformMtrx[0][2]);
		XMLMatrix->SetDoubleAttribute ("z2", Helper.nativeTransformMtrx[1][2]);
		XMLMatrix->SetDoubleAttribute ("z3", Helper.nativeTransformMtrx[2][2]);
		Root->AddNestedElement(XMLMatrix);
	}
	if(anOrientation)
	{
		vtkXMLDataElement *XMLAnOr = vtkXMLDataElement::New ();
		XMLAnOr->SetName("anamori");
		XMLAnOr->SetCharacterData (Helper.anatomicalOrientation, strlen(Helper.anatomicalOrientation)+1);
		Root->AddNestedElement(XMLAnOr);
	}
	
	vtkXMLDataElement *Slices=vtkXMLDataElement::New ();
	Slices->SetIntAttribute ("orientation", Helper.xyzOrientation);
	Slices->SetIntAttribute ("offset", Helper.offset);
	Slices->SetIntAttribute ("timestamp",(int) vtkTimerLog::GetUniversalTime());
	Slices->SetName("slicesprefix");
	Slices->SetCharacterData (Helper.prefix, strlen(Helper.prefix)+1);
	Root->AddNestedElement(Slices);

	std::stringstream FileName;
	FileName<<Helper.directory<<"volume.xml";
	vtkXMLUtilities::WriteElementToFile (Root, FileName.str().c_str());
}
