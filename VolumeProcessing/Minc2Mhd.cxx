/*=========================================================================

Program:   Minc2Mhd : converts .mnc images to mhd image format
Module:    vtkSurface
Language:  C++
Date:      2010/04
Auteur:   Sebastien Valette

=========================================================================*/
// .NAME Minc2Mhd
// .SECTION Description

#include <sstream>
#include <vtkImageData.h>
#include <vtkMetaImageWriter.h>
#include <vtkMINCImageReader.h>

int main( int argc, char *argv[] )
{

	if (argc<3)
	{
		cout<<"Usage : Minc2Mhd file.mnc file.mhd"<<endl;
		exit(1);
	}

	// Load Volume
	cout <<"load : "<<argv[1]<<endl;

	vtkMINCImageReader *Reader=vtkMINCImageReader::New();
	Reader->SetFileName(argv[1]);
	Reader->Update();
	vtkMetaImageWriter *Writer=vtkMetaImageWriter::New();
	Writer->SetInputData(Reader->GetOutput());
	Writer->SetFileName(argv[2]);
	Writer->Write();
}
