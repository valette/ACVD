/*
 * Authors: R. Kechichian
 *          S. Valette
 */

#include <string>
#include <vector>
#include <iomanip>
#include <iostream>
#include <algorithm>

using std::string ;
using std::vector ;
using std::ostream ;
using std::cout ;
using std::cerr ;
using std::setw ;

#include <vtkDirectory.h>
#include <vtkImageData.h>
#include <vtkSmartPointer.h>
#include <vtkImageReader2.h>
#include <vtkMetaImageReader.h>
#include <vtkMetaImageWriter.h>
#include <vtkImageReader2Factory.h>


void printSyntax(ostream& os)
{
	os << "img2raw [-c] <image> ... <directory> ... -o <output filename>\n" ;
	os << "Creates a compressed or uncompressed Meta Image (MHD + RAW)\n" ;
	os << "from input images and/or directories of images specified in any order.\n" ;
	os << "\t* Supported image formats:\n" ;
	os << "\t\tPNG, PNM, TIFF, BMP, SLC, JPEG, GE Signa, MINC, Meta Image (MHD + RAW)\n" ;
	os << "\t* Ouput volume is of dimensions and voxel data type as to accommodate all input.\n" ;
	os << "\t* Z-order of output volume is determined by input order.\n" ;
	os << "\t* Doesn't recurse subdirectories.\n" ;
}

int main(int argc, char *argv[])
{
	// Parse and validate input
	vector<string> arguments ;
	for (int i = 0 ; i < argc ; ++i )
		arguments.push_back(argv[i]) ;

	vector<string>::iterator itOptO = find(arguments.begin(), arguments.end(), "-o") ;
	if (itOptO == arguments.end()) {
		cerr << "Required option -o missing.\n" ;
		printSyntax(cerr) ;
		return 1 ;
	}

	vector<string>::iterator itArgO = itOptO ;
	++itArgO ;
	if (itArgO == arguments.end()) {
		cerr << "Required output raw image filename missing.\n" ;
		printSyntax(cerr) ;
		return 2 ;
	}

	string outputFileName = *itArgO ;
	string mhdFileName, rawFileName ;
	int pointPos = outputFileName.find_last_of('.') ;
	if (pointPos < outputFileName.length()) {
		string outputFileExtension = outputFileName.substr(pointPos + 1, outputFileName.length()) ;
		if ((outputFileExtension == "raw") || (outputFileExtension == "mhd")) {
			string fileName = outputFileName.substr(0, pointPos) ;
			rawFileName = fileName + ".raw" ;
			mhdFileName = fileName + ".mhd" ;
		}
		else {
			rawFileName = outputFileName + ".raw" ;
			mhdFileName = outputFileName + ".mhd" ;
		}
	}
	else {
		rawFileName = outputFileName + ".raw" ;
		mhdFileName = outputFileName + ".mhd" ;
	}

	bool compressRAW = false ;
	vector<string>::iterator itOptC = find(arguments.begin(), arguments.end(), "-c") ;
	if (itOptC != arguments.end())
		compressRAW = true ;

	arguments.erase(itArgO) ;
	itOptO = find(arguments.begin(), arguments.end(), "-o") ;
	if (itOptO != arguments.end())
		arguments.erase(itOptO) ;
	itOptC = find(arguments.begin(), arguments.end(), "-c") ;
	if (itOptC != arguments.end())
		arguments.erase(itOptC) ;
	vector<string>::iterator itProg = find(arguments.begin(), arguments.end(), argv[0]) ;
	arguments.erase(itProg) ;

	if (!arguments.size()) {
		cerr << "No input directories or images to process.\n" ;
		printSyntax(cerr) ;
		return 3 ;
	}

	cout << "Reading input... " ;
	cout.flush() ;

	// Create image reader factory and register Meta Image reader with it
	vtkSmartPointer<vtkImageReader2Factory> imageReaderFactory = vtkSmartPointer<vtkImageReader2Factory>::New() ;
	vtkSmartPointer<vtkMetaImageReader> metaImageReader = vtkSmartPointer<vtkMetaImageReader>::New() ;
	imageReaderFactory->RegisterReader(metaImageReader) ;

	int volumeWidth = 0, volumeHeight = 0, volumeDepth = 0 ;
	int volumeScalarSize = 0, volumeScalarComponents = 0 ;

	// Read input directories and/or images
	vector<vtkImageReader2*> imageReaders ;
	for (int i = 0 ; i < arguments.size() ; ++i) {
		vtkSmartPointer<vtkDirectory> directory = vtkSmartPointer<vtkDirectory>::New() ;
		if (directory->FileIsDirectory(arguments[i].c_str())) {
			if (!directory->Open(arguments[i].c_str())) {
				cerr << "Cannot open directory " << arguments[i] << "; skipping it.\n" ;
				continue ;
			}
			int numberOfFiles = directory->GetNumberOfFiles() ;
			vector<string> fileNames ;
			for (int j = 0 ; j < numberOfFiles ; ++j) {
				if ((string(".") == directory->GetFile(j)) || (string("..") == directory->GetFile(j)))
					continue ;
				string fileName = arguments[i] + string("/") + directory->GetFile(j) ;
				vtkSmartPointer<vtkDirectory> auxDirectory = vtkSmartPointer<vtkDirectory>::New() ;
				if (auxDirectory->FileIsDirectory(fileName.c_str())) {
					cerr << "Found sub-directory " << fileName << "; skipping it.\n" ;
					continue ;
				}
				fileNames.push_back(fileName) ;
			}
			sort(fileNames.begin(), fileNames.end()) ;
			for (int j = 0 ; j < fileNames.size() ; ++j) {
				vtkImageReader2* imageReader = imageReaderFactory->CreateImageReader2(fileNames[j].c_str()) ;
				if (!imageReader) {
					cerr << "Cannot load file " << fileNames[j] << " as an image file; skipping it.\n" ;
					continue ;
				}
				else {
					imageReader->SetFileName(fileNames[j].c_str()) ;
					imageReader->Update() ;
					int imageDimensions[3] ;
					imageReader->GetOutput()->GetDimensions(imageDimensions) ;
					if (imageDimensions[0] > volumeWidth)
						volumeWidth = imageDimensions[0] ;
					if (imageDimensions[1] > volumeHeight)
						volumeHeight = imageDimensions[1]  ;
					volumeDepth += imageDimensions[2] ;
					int imageScalarSize = imageReader->GetOutput()->GetScalarSize() ;
					if (imageScalarSize > volumeScalarSize)
						volumeScalarSize = imageScalarSize ;
					int imageScalarComponents = imageReader->GetOutput()->GetNumberOfScalarComponents() ;
					if (imageScalarComponents > volumeScalarComponents)
						volumeScalarComponents = imageScalarComponents ;
					imageReaders.push_back(imageReader) ;
				}
			}
		}
		else {
			vtkImageReader2* imageReader = imageReaderFactory->CreateImageReader2(arguments[i].c_str()) ;
			if (!imageReader) {
				cerr << "Cannot load file " << arguments[i] << " as an image file; skipping it.\n" ;
				continue ;
			}
			else {
				imageReader->SetFileName(arguments[i].c_str()) ;
				imageReader->Update() ;
				int imageDimensions[3] ;
				imageReader->GetOutput()->GetDimensions(imageDimensions) ;
				if (imageDimensions[0] > volumeWidth)
					volumeWidth = imageDimensions[0] ;
				if (imageDimensions[1] > volumeHeight)
					volumeHeight = imageDimensions[1] ;
				volumeDepth += imageDimensions[2] ;
				int imageScalarSize = imageReader->GetOutput()->GetScalarSize() ;
				if (imageScalarSize > volumeScalarSize)
					volumeScalarSize = imageScalarSize ;
				int imageScalarComponents = imageReader->GetOutput()->GetNumberOfScalarComponents() ;
				if (imageScalarComponents > volumeScalarComponents)
					volumeScalarComponents = imageScalarComponents ;
				imageReaders.push_back(imageReader) ;
			}
		}
		cout << '\r' << "Reading input... " << setw(3) << static_cast<int>(100.0 * (static_cast<double>(i+1) / arguments.size())) << "% " ;
		cout.flush() ;
	}

	cout << "done.\n" ;
	cout.flush() ;

	if (!imageReaders.size()) {
		cerr << "No valid input has been read, terminating.\n" ;
		return 4 ;
	}

	// Create and construct volume
	vtkSmartPointer<vtkImageData> volumeImageData = vtkSmartPointer<vtkImageData>::New() ;
	volumeImageData->SetDimensions(volumeWidth, volumeHeight, volumeDepth) ;
	volumeImageData->SetExtent(0, volumeWidth-1, 0, volumeHeight-1, 0, volumeDepth-1) ;

	int type;
	if (volumeScalarSize == 1)
		type = VTK_UNSIGNED_CHAR;
	else if (volumeScalarSize == 2)
		type = VTK_UNSIGNED_SHORT;
	else if (volumeScalarSize == 4)
		type = VTK_UNSIGNED_INT;
	else
		type = VTK_FLOAT;
	
	volumeImageData->AllocateScalars(type, volumeScalarComponents);

	cout << "Constructing volume... ";
	cout.flush() ;

	int volumeX = 0, volumeY = 0, volumeZ = 0 ;
	char *iP = (char *) volumeImageData->GetScalarPointer(), *iP2;
	int dimensions[3];
	volumeImageData->GetDimensions(dimensions);

	int size = volumeImageData->GetScalarSize() 
		* dimensions[0] * dimensions[1]
		* volumeImageData->GetNumberOfScalarComponents();

	for (int i = 0 ; i < imageReaders.size() ; ++i) {
		vtkImageData* imageData = imageReaders[i]->GetOutput();
		cout << '\r' << "Constructing volume... " << setw(3) << static_cast<int>(100.0 * (static_cast<double>(i) / imageReaders.size())) << "% " ;
		cout.flush() ;
		iP2 = (char *) imageData->GetScalarPointer();
		for (int j = 0; j < size ; j++) {
			*(iP++) = *(iP2++);
		}
		imageReaders[i]->Delete() ;
	}

	cout << "done.\n" ;
	cout.flush() ;

	cout << "Writing volume data in " << rawFileName << "... " ;
	cout.flush() ;

	// Write out volume
	vtkSmartPointer<vtkMetaImageWriter> volumeWriter = vtkSmartPointer<vtkMetaImageWriter>::New() ;
	volumeWriter->SetFileName(mhdFileName.c_str()) ;
	volumeWriter->SetRAWFileName(rawFileName.c_str()) ;
	volumeWriter->SetInputData(volumeImageData) ;
	volumeWriter->SetCompression(compressRAW) ;
	volumeWriter->Write() ;

	cout << "done.\n" ;
	cout.flush() ;

	cout << "Output image: " << mhdFileName << '\n' ;
	cout << "\tVolume dimensions: width = " << volumeWidth << ", height = " << volumeHeight << ", depth = " << volumeDepth << '\n' ;
	cout << "\tVoxels have " << volumeScalarComponents << " components of size " << volumeScalarSize << " bytes.\n" ;
}
