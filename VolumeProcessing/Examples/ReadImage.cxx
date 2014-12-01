/*=========================================================================

Program:   tire : folds an image to a tire-like 3D model
Module:    vtkSurface
Language:  C++
Date:      2010/03
Auteur:   Sebastien Valette

=========================================================================*/
// .NAME meshviewer 
// .SECTION Description

#include <vtkPNGReader.h>
#include <vtkImageViewer2.h>
#include <vtkImageReader.h>
#include <vtkImageData.h>
#include <vtkRenderWindowInteractor.h>

int main( int argc, char *argv[] )
{
	vtkImageData *Image=0;

	// Load and Display Image
	cout <<"load : "<<argv[1]<<endl;

	vtkPNGReader *Reader=vtkPNGReader::New();
	Reader->SetFileName(argv[1]);
	Reader->Update();
	Image=Reader->GetOutput();

	/*vtkImageViewer2 *Viewer=vtkImageViewer2::New();
	Viewer->SetInput(Image);

	vtkRenderWindowInteractor *Interactor=vtkRenderWindowInteractor::New();
	Viewer->SetupInteractor(Interactor);
	Interactor->Start();
	Viewer->Render();
*/
	cout<<"Scalar Type : "<<Image->GetScalarType()<<endl;
	cout<<"Valeur "<<*(unsigned short*)Image->GetScalarPointer(220,272,0)<<endl;;
}
