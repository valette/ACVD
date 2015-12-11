/*=========================================================================

Program:   VolumeSubsample :subsamples a volume
Module:    vtkSurface
Language:  C++
Date:      2010/04
Auteur:   Sebastien Valette

=========================================================================*/
// .NAME VolumeSubsample2Mhd
// .SECTION Description

#include <sstream>
#include <vtkImageData.h>
#include <vtkMetaImageWriter.h>
#include <vtkMetaImageReader.h>
#include <vtkImageResample.h>

int main( int argc, char *argv[] )
{

	if (argc<3)
	{
		cout<<"Usage : Minc2Mhd in.mhd out.mhd factor"<<endl;
		exit(1);
	}

	// Load Volume
	cout <<"load : "<<argv[1]<<endl;

	vtkMetaImageReader *Reader=vtkMetaImageReader::New();
	Reader->SetFileName(argv[1]);
	Reader->Update();

	vtkImageResample *Resampler=vtkImageResample::New();
	Resampler->SetInputConnection(Reader->GetOutputPort());
	Resampler->SetInterpolationModeToNearestNeighbor();
	for (int i=0;i<3;i++)
		Resampler->SetAxisMagnificationFactor(i, atof(argv[3]));
	cout<<"Resampling...";
	Resampler->Update();
	vtkImageData *Image = Resampler->GetOutput();
	cout<<"done."<<endl;
	int Dimensions[3];
	Image->GetDimensions(Dimensions);
	cout<<"New image dimensions : "<<Dimensions[0]<<" "<<Dimensions[1]<<" "<<Dimensions[2]<<endl;

	vtkMetaImageWriter *Writer = vtkMetaImageWriter::New();
	Writer->SetInputData(Image);
	Writer->SetFileName(argv[2]);
	Writer->Write();
}
