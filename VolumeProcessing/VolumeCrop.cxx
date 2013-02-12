/*=========================================================================

Program:   VolumeCrop :crops a volume
Module:    vtkSurface
Language:  C++
Date:      2012/09
Auteur:   Sebastien Valette

=========================================================================*/
// .NAME VolumeCrop
// .SECTION Description

#include <sstream>
#include <vtkImageData.h>
#include <vtkMetaImageWriter.h>
#include <vtkMetaImageReader.h>
#include <vtkExtractVOI.h>

int main( int argc, char *argv[] )
{

	if (argc<9)
	{
		cout<<"Usage : VolumeCrop in.mhd out.mhd xmin xmax ymin ymax zmin zmax"<<endl;
		exit(1);
	}

	// Load Volume
	cout <<"load : "<<argv[1]<<endl;


	// read image
	vtkMetaImageReader *Reader=vtkMetaImageReader::New();
	Reader->SetFileName(argv[1]);
	Reader->Update();

	// extract VOI
	int VOI[6];
	for (int i=0;i<6;i++)
	{
		VOI[i]=atoi(argv[i+3]);
	}
	vtkExtractVOI *Extract=vtkExtractVOI::New();
	Extract->SetVOI(VOI);
	Extract->SetInput(Reader->GetOutput());
	Extract->Update();

	//save VOI
	vtkImageData *Image=Extract->GetOutput();
	cout<<"done."<<endl;
	int Dimensions[3];
	Image->GetDimensions(Dimensions);
	cout<<"New image dimensions : "<<Dimensions[0]<<" "<<Dimensions[1]<<" "<<Dimensions[2]<<endl;

	vtkMetaImageWriter *Writer=vtkMetaImageWriter::New();
	Writer->SetInput(Image);
	Writer->SetFileName(argv[2]);
	Writer->Write();
}
