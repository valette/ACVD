/*=========================================================================

Program:   VolumeMedian : processes median filtering on input volume
Module:    vtkSurface
Language:  C++
Date:      2011/09
Auteur:   Sebastien Valette

=========================================================================*/
// .NAME VolumeMedian
// .SECTION Description

#include <vtkImageData.h>
#include <vtkMetaImageReader.h>
#include <vtkMetaImageWriter.h>
#include <vtkImageMedian3D.h>

int main( int argc, char *argv[] )
{
	if (argc<3)
	{
		cout<<"Usage : VolumeMedian file.mhd radius"<<endl;
		exit(1);
	}

	// Load Volume
	cout <<"load : "<<argv[1]<<endl;

	vtkMetaImageReader *Reader=vtkMetaImageReader::New();
	Reader->SetFileName(argv[1]);
	Reader->Update();
	vtkImageMedian3D *Median=vtkImageMedian3D::New();
	Median->SetInputConnection(Reader->GetOutputPort());
	int Radius=atoi(argv[2]);
	Median->SetKernelSize(Radius,Radius,Radius);
	Median->Update();
	vtkMetaImageWriter *Writer=vtkMetaImageWriter::New();
	Writer->SetInputData(Median->GetOutput());
	Writer->SetFileName("output.mhd");
	Writer->Write();
}
