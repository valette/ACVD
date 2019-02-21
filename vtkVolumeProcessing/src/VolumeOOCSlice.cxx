/*=========================================================================

Program:   VolueSection : extracts a .png slice from a volume in an
			out-of-core manner
Module:    vtkSurface
Language:  C++
Date:      2015/03
Auteur:   Sebastien Valette / J. Eckert

=========================================================================*/
// .NAME VolueSection
// .SECTION Description

#include "vtkOOCMetaImageReader.h"
#include <vtkImageData.h>
#include <vtkJPEGWriter.h>
#include <vtkPNGWriter.h>
#include <vtkPointData.h>
#include <vtkImageShiftScale.h>
#include <vtkSmartPointer.h>
#include <vtkImageFlip.h>

int main( int argc, char *argv[] )
{
	if (argc < 3) {
		cout<<"Usage : volumeOOCSlice file slice"<<endl;
		exit(1);
	}

	vtkOOCMetaImageReader *Reader = vtkOOCMetaImageReader::New();

	Reader->SetFileName(argv[1]);

	int z = atoi(argv[2]);

	Reader->SetZMin(z);
	Reader->SetZMax(z);
	Reader->Update();

	vtkImageFlip *flip = vtkImageFlip::New();
	flip->SetInput(Reader->GetOutput());
	flip->SetFilteredAxis(1);
	flip->Update();
	vtkImageData *Image = flip->GetOutput();

	int scalarType = Image->GetScalarType();
	vtkPNGWriter* writer = vtkPNGWriter::New();

	int NumComponents;
	switch (Image->GetScalarType()) {
	case 2:
	case 15:
	case 3:
		NumComponents = 1;
		break;
	case 4:
	case 5:
		NumComponents = 2;
		break;
	default:
		NumComponents = 4;
		break;
	}

	int Dimensions[3];
	Image->GetDimensions(Dimensions);
	Dimensions[2] = 1;

	vtkImageData *PNGSlice = vtkImageData::New();
	PNGSlice->SetNumberOfScalarComponents(NumComponents);
	PNGSlice->SetScalarTypeToUnsignedChar();
	PNGSlice->SetDimensions(Dimensions);
	PNGSlice->AllocateScalars();		
	unsigned char* PNGPointer = (unsigned char *) PNGSlice->GetScalarPointer();

	if (Image->GetNumberOfScalarComponents() > 1) {
		writer->SetInput(Image);
	} else {
		writer->SetInput(PNGSlice);

		vtkImageData *Output;
		unsigned char* SlicePointer;
		int numPixels;

		switch (Image->GetScalarType()) {
		case 2:
		case 15:
		case 3:
		case 4:
		case 5:
			Output = Image;
			break;
		default:
			vtkImageShiftScale *Cast = vtkImageShiftScale::New();
			Cast->SetOutputScalarTypeToFloat ();
			Cast->SetInput(Image);
			Cast->Update();
			Output = Cast->GetOutput();
			break;
		}

		SlicePointer = (unsigned char *) Output->GetScalarPointer();
		numPixels = Dimensions[0] * Dimensions[1] * NumComponents;
		for (int j = 0; j != numPixels; j++) {
			PNGPointer[j] = SlicePointer[j];
		}
	}

	writer->SetFileName("slice.png");
	writer->Write();

	double range[2];
	Image->GetScalarRange(range);
    std::ofstream out("range.txt", std::ios::out | std::ofstream::trunc);
    out << range[0] << " " << range[1];
    out.close();

}
