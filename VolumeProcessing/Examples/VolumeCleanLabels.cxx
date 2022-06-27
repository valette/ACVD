/*=========================================================================

Program:   VolumeCleanLabels : keeps only biggest labeled connected components
Module:    vtkSurface
Language:  C++
Date:      2012/04
Auteur:   Sebastien Valette

=========================================================================*/
// .NAME VolumeCleanLabels
// .SECTION Description

#include <vtkMetaImageReader.h>
#include <vtkMetaImageWriter.h>
#include "vtkImageDataCleanLabels.h"

int main( int argc, char *argv[] )
{
	if (argc<2)
	{
		cout<<"Usage : VolumeCleanLabels file.mhd"<<endl;
		exit(1);
	}

	// Load Volume
	cout <<"load : "<<argv[1]<<endl;

	vtkMetaImageReader *Reader=vtkMetaImageReader::New();
	Reader->SetFileName(argv[1]);
	Reader->Update();
	vtkImageDataCleanLabels *Cleaner=vtkImageDataCleanLabels::New();
	Cleaner->SetInputConnection(Reader->GetOutputPort());
	Cleaner->Update();
	vtkMetaImageWriter *Writer=vtkMetaImageWriter::New();
	Writer->SetInputConnection(Cleaner->GetOutputPort());
	Writer->SetFileName("output.mhd");
	Writer->Write();
}
