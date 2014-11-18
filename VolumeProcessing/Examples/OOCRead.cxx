/*=========================================================================

Program:   png2raw : loads png images and generates a .raw volume
Module:    vtkSurface
Language:  C++
Date:      2010/05
Auteur:   Sebastien Valette

=========================================================================*/
// .NAME png2raw 
// .SECTION Description

#include "vtkOOCMetaImageReader.h"
#include "ResampledVolumeGenerator.h"
#include "VolumeRenderWindow.h"

int main( int argc, char *argv[] )
{
	if (argc<8)
	{
		cout<<"Usage : png2raw file x1 x2 y1 y2 z1 z2"<<endl;
		exit(1);
	}

	vtkOOCMetaImageReader *Reader=vtkOOCMetaImageReader::New();

	Reader->SetFileName (argv[1]);

	int VOI[6];

	for (int i=0;i<6;i++)
		VOI[i]=atoi(argv[i+2]);

	Reader->SetDataVOI (VOI);
	Reader->Update();
	
	VolumeRenderWindow *Window=VolumeRenderWindow::New();
	vtkImageData *Image=Reader->GetOutput();
	
	int Dimensions[3];
	Image->GetDimensions(Dimensions);
	
	cout<<"Image dimensions : ";
	for (int i=0;i<3;i++)
		cout<<Dimensions[i]<<" ";
	cout<<endl;
	
	cout<<"Scalar Type : "<<Image->GetScalarType()<<endl;
	cout<<"Number of Components : "<<Image->GetNumberOfScalarComponents()<<endl;
	
	Window->SetInput(Image);
	Window->Render();
	Window->Interact();
	/*
	
	ResampledVolumeGenerator *rvg = new ResampledVolumeGenerator();
	rvg->SetFileName(argv[1]);
	rvg->SetDivisionFactors(5,5,5);
	rvg->SetResamplingFactors(0.1,0.1,0.1);
	rvg->UnthreadedExecute();*/
}
