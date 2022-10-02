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
#include <vtkExtractVOI.h>
#include <vtkNIFTIImageWriter.h>

#include "vtkRobustImageReader.h"

int main( int argc, char *argv[] )
{

	if (argc<9)
	{
		cout<<"Usage : VolumeCrop in.mhd out.nii.gz xmin xmax ymin ymax zmin zmax"<<endl;
		exit(1);
	}

	// Load Volume
	cout <<"load : "<<argv[1]<<endl;

	// read image
	vtkRobustImageReader *Reader=vtkRobustImageReader::New();
	Reader->SetFileName(argv[1]);
	Reader->Update();

	// extract VOI
	int VOI[6];
	for (int i=0;i<6;i++)
		VOI[i]=atoi(argv[i+3]);

	vtkExtractVOI *Extract=vtkExtractVOI::New();
	Extract->SetVOI(VOI);
	Extract->SetInputData(Reader->GetOutput());
	Extract->Update();

	//save VOI
	vtkImageData *Image=Extract->GetOutput();
	cout<<"done."<<endl;
	int Dimensions[3];
	Image->GetDimensions(Dimensions);
	cout<<"New image dimensions : "<<Dimensions[0]<<" "<<Dimensions[1]<<" "<<Dimensions[2]<<endl;

	double origin[ 3 ], spacing[ 3 ];
	Reader->GetOutput()->GetOrigin( origin );
	Reader->GetOutput()->GetSpacing( spacing );

	for ( int i = 0; i < 3; i++ )
		origin[ i ] = origin[ i ] + spacing[ i ] * VOI[ 2 * i ];

	Image->SetOrigin( origin );
	vtkNIFTIImageWriter *Writer=vtkNIFTIImageWriter::New();
	Writer->SetInputData(Image);
	Writer->SetFileName(argv[2]);
	Writer->Write();

}
