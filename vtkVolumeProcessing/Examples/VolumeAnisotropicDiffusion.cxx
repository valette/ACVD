/*=========================================================================

Program:   VolumeAnisotropicDiffusion : processes median filtering on input volume
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
#include <vtkImageAnisotropicDiffusion3D.h>

int main( int argc, char *argv[] )
{
	if (argc<3)
	{
		cout<<"Usage : VolumeAnisotropicDiffusion file.mhd numberofiterations []"<<endl;
		exit(1);
	}

	// Load Volume
	cout <<"load : "<<argv[1]<<endl;

	vtkMetaImageReader *Reader=vtkMetaImageReader::New();
	Reader->SetFileName(argv[1]);
	Reader->Update();
	vtkImageAnisotropicDiffusion3D *Diffusion=vtkImageAnisotropicDiffusion3D::New();
	Diffusion->SetInputData(Reader->GetOutput());
	int NumberOfIterations=atoi(argv[2]);
	Diffusion->SetNumberOfIterations(NumberOfIterations);

	// Parse optionnal arguments
	int ArgumentsIndex=3;
	while (ArgumentsIndex<argc)
	{
		if (strcmp(argv[ArgumentsIndex],"-t")==0)
		{
			Diffusion->SetDiffusionThreshold(atof(argv[ArgumentsIndex+1]));
			cout<<"Diffusion threshold="<<atof(argv[ArgumentsIndex+1])<<endl;
		}

		if (strcmp(argv[ArgumentsIndex],"-f")==0)
		{
			Diffusion->SetDiffusionFactor(atof(argv[ArgumentsIndex+1]));
			cout<<"Diffusion factor="<<atof(argv[ArgumentsIndex+1])<<endl;
		}
		ArgumentsIndex+=2;
	}

	Diffusion->Update();
	vtkMetaImageWriter *Writer=vtkMetaImageWriter::New();
	Writer->SetInputData(Diffusion->GetOutput());
	Writer->SetFileName("output.mhd");
	Writer->Write();
}
