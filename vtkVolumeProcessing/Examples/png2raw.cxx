/*=========================================================================

Program:   png2raw : loads png images and generates a .raw volume
Module:    vtkSurface
Language:  C++
Date:      2010/05
Auteur:   Sebastien Valette

=========================================================================*/
// .NAME png2raw
// .SECTION Description

#include <vtkPNGReader.h>
#include <vtkImageViewer2.h>
#include <vtkImageReader.h>
#include <vtkImageData.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkMetaImageWriter.h>
#include <sstream>

int main( int argc, char *argv[] )
{
	if (argc<4)
	{
		cout<<"Usage : png2raw fileprefix initial_index number_of_images"<<endl;
		exit(1);
	}

	int StartImage=atoi(argv[2]);
	int NumberOfImages=atoi(argv[3]);

	fstream  binary_file("volume.raw",ios::out|ios::binary);
	int CurrentImage=StartImage;

	for (int i=0;i<NumberOfImages;i++)
	{
		std::stringstream Name;
		Name<<argv[1]<<CurrentImage<<".png";

		// Load and Display Image
		cout <<"load : "<<Name.str()<<endl;

		vtkPNGReader *Reader=vtkPNGReader::New();
		Reader->SetFileName(Name.str().c_str());
		Reader->Update();
		vtkImageData *Image=Reader->GetOutput();

		cout<<"Scalar Type : "<<Image->GetScalarType()<<endl;
//		cout<<"Valeur "<<*(unsigned short*)Image->GetScalarPointer(220,272,0)<<endl;;
		char *Pointer=(char *) Image->GetScalarPointer(0,0,0);

		int Dimensions[3];
		Image->GetDimensions(Dimensions);
		binary_file.write(Pointer, sizeof(unsigned char)*Dimensions[0]*Dimensions[1]*Image->GetNumberOfScalarComponents());
		Reader->Delete();
		CurrentImage++;
	}
	binary_file.close();
	/*
	vtkMetaImageWriter *Writer=vtkMetaImageWriter::New();
	Writer->SetCompression(0);
	Writer->SetInput(Image);
	Writer->SetFileName("output.mhd");
	Writer->Write();

	vtkImageViewer2 *Viewer=vtkImageViewer2::New();
	Viewer->SetInput(Image);

	vtkRenderWindowInteractor *Interactor=vtkRenderWindowInteractor::New();
	Viewer->SetupInteractor(Interactor);
	Interactor->Start();
	Viewer->Render();
	*/
}
