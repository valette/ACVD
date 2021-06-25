#ifndef __vtkRobustImageReader_h
#define __vtkRobustImageReader_h

#include <sstream>

#include <vtkImageData.h>
#include <vtkImageFlip.h>
#include <vtkImageReader2Factory.h>
#include <vtkImageShiftScale.h>
#include <vtkMatrix4x4.h>
#include <vtkMetaImageReader.h>
#include <vtkNIFTIImageReader.h>
#include <vtkObjectFactory.h>

// v2 : shift and scale values read by nifti reader when needed
// v1

class vtkRobustImageReader : public vtkObject 
{

public :

	static vtkRobustImageReader *New() {
		vtkObject* ret = vtkObjectFactory::CreateInstance("vtkRobustImageReader");
		if(ret) {
			return (vtkRobustImageReader*)ret;
		}

		return (new vtkRobustImageReader);
	}

	vtkTypeMacro(vtkRobustImageReader,vtkObject);

	void Update() {
		vtkImageReader2Factory *imageReaderFactory = vtkImageReader2Factory::New();
		vtkMetaImageReader *metaImageReader = vtkMetaImageReader::New();
		imageReaderFactory->RegisterReader(metaImageReader);
		vtkNIFTIImageReader *niftiiImageReader = vtkNIFTIImageReader::New();
		imageReaderFactory->RegisterReader(niftiiImageReader);

		// Create a reader for image and try to load it
		vtkImageReader2* Reader = imageReaderFactory->CreateImageReader2(FileName) ;

		if (!Reader) {
			std::cerr << "Cannot load file " << FileName << " as an image file; terminating.\n" ;
			exit(5) ;
		}

		Reader->SetFileName(FileName);
		Reader->Update();
		Output = Reader->GetOutput();

		bool flip[3] = {false, false, false};
		double Origin[3];
		double Spacing[3];
		int Dimensions[3];

		if (strcmp (Reader->GetClassName(), "vtkMetaImageReader") == 0) {
			bool matrixFound = false;
			bool anOrientation = false;
			std::string matrixPrefix;
			std::string anatomicalOrientation;

			std::ifstream is(FileName);
			std::string line, prefix;
			double transformationMatrix[9] ;

			while(is) {
				getline(is, line);
				std::istringstream iss(line);
				iss >> prefix;
				if ((prefix=="TransformMatrix") || (prefix == "Orientation") || (prefix == "Rotation")) {
					matrixPrefix = prefix;
					iss >> prefix; // eat the "=" sign
					for (int i = 0; i < 9; i++) {
						iss >> transformationMatrix[i];
					}

					for (int i = 0; i < 3; i++) {
						if (transformationMatrix[4 * i] < 0) {
							flip[i] = true;
						}
					}
				}
			}
		}

		if (strcmp (Reader->GetClassName(), "vtkNIFTIImageReader") == 0) {
		    vtkNIFTIImageReader *niftiReader = (vtkNIFTIImageReader *) Reader;

            vtkMatrix4x4 *qForm = niftiReader->GetQFormMatrix();
			if (!qForm) {
                qForm = vtkMatrix4x4::New();
            }

            for (int i = 0; i < 3; i++) {
				if (qForm->GetElement(i, i) < 0) {
					flip[i] = true;
				}
				Origin[i] = qForm->GetElement(i,3);
			}


            if (niftiReader->GetQFac() < 0) {
				for (int i = 0; i < 2; i++) {
					flip[i] = !flip[i];
					Origin[i] = -Origin[i];
				}
            }

			Output->SetOrigin(Origin);
		}

		if (!flip[0] && !flip[1] && !flip[2]) {
			return;
		}

		Output->GetSpacing(Spacing);
		Output->GetDimensions(Dimensions);

        for (int i = 0; i < 3; i++) {
            if (!flip[i]) continue;
            std::cout << "Warning : RobustReader flipping dimension " << i << std::endl;
            vtkImageFlip *flip = vtkImageFlip::New();
            flip->SetInputData(Output);
            flip->SetFilteredAxis (i);
            flip->Update();
            Output = flip->GetOutput();
			Output->GetOrigin(Origin);
            Origin[i] = Origin[i] - Spacing[i] * ( Dimensions[i] - 1);
			Output->SetOrigin(Origin);
        }

		if (strcmp (Reader->GetClassName(), "vtkNIFTIImageReader") == 0) {
		    vtkNIFTIImageReader *niftiReader = (vtkNIFTIImageReader *) Reader;

			double slope = niftiReader->GetRescaleSlope();
			double intercept = niftiReader->GetRescaleIntercept();

			if ( ( slope != 1.0 ) || ( intercept != 0 ) ) {

				vtkImageShiftScale *shiftScale = vtkImageShiftScale::New();
				shiftScale->SetShift( intercept );
				shiftScale->SetScale( slope );
				shiftScale->SetInputData( Output );
				shiftScale->Update();
				Output = shiftScale->GetOutput();
				std::cout << "Warning : RobustReader shifting nifti volume values" << std::endl;

			}

		}

	}

	vtkGetObjectMacro(Output, vtkImageData)

	vtkGetMacro(FileName, char*)
	vtkSetMacro(FileName, char*)

protected :
	vtkImageData *Output;

	char *FileName;

	vtkRobustImageReader() {
		Output = 0;
		FileName = 0;
	};
};

#endif

