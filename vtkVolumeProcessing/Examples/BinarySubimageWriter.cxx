#include "BinarySubimageWriter.h"

//ENCORE A FAIRE: Bien définir les destructeurs de mes classes afin d'y inclure les 'deletes' de chaque élement
BinarySubimageWriter::BinarySubimageWriter(){
	_image = vtkImageData::New();
	_subimageExtent = new int[6];
	_imageDimensions = new int[3];
	_pixelSize = 0;
	//_Filename = "";
}

BinarySubimageWriter::~BinarySubimageWriter(){
	_image->Delete();
	delete _Filename;
	delete _subimageExtent;
	delete _imageDimensions;
}

//static BinarySubimageWriter* GetInstance(){
//	if(instance == NULL){
//		instance = new BinarySubimageWriter();
//	}
//	return instance;
//}

void BinarySubimageWriter::SetNumberOfChannels ( int numberOfChannels ) {
	_numberOfChannels = numberOfChannels;
}

void BinarySubimageWriter::SetPixelSizeAndType ( int pixelType ){
	switch(pixelType)
	{
	case VTK_SIGNED_CHAR:
		_pixelType = VTK_SIGNED_CHAR;
		_pixelSize = sizeof(char);
		break;
	case VTK_UNSIGNED_CHAR:
		_pixelType = VTK_UNSIGNED_CHAR;
		_pixelSize = sizeof(unsigned char);
		break;
	case VTK_SHORT:
		_pixelType = VTK_SHORT;
		_pixelSize = sizeof(short);
		break;
	case VTK_UNSIGNED_SHORT:
		_pixelType = VTK_UNSIGNED_SHORT;
		_pixelSize = sizeof(unsigned short);
		break;
	case VTK_INT:
		_pixelType = VTK_INT;
		_pixelSize = sizeof(int);
		break;
	case VTK_UNSIGNED_INT:
		_pixelType = VTK_UNSIGNED_INT;
		_pixelSize = sizeof(unsigned int);
		break;
	case VTK_LONG:
		_pixelType = VTK_LONG;
		_pixelSize = sizeof(long);
		break;
	case VTK_UNSIGNED_LONG:
		_pixelType = VTK_UNSIGNED_LONG;
		_pixelSize = sizeof(unsigned long);
		break;
	case VTK_FLOAT:
		_pixelType = VTK_FLOAT;
		_pixelSize = sizeof(float);
		break;
	case VTK_DOUBLE:
		_pixelType = VTK_DOUBLE;
		_pixelSize = sizeof(double);
		break;
	}
}

void BinarySubimageWriter::AssignFilename ( char* originalFilename ){
	//int i;
	//long l = static_cast<long>( strlen( originalFilename ) );
	//_Filename = new char[l+4];

	//for(i=l; i>=0; i--){
	//	if( originalFilename[i] == '.' ){
	//		strcpy(_Filename, originalFilename);
	//		_Filename[i] = '_';
	//		_Filename[i+1] = 'a';
	//		_Filename[i+2] = 'l';
	//		_Filename[i+3] = 't';
	//		_Filename[i+4] = '.';
	//		_Filename[i+5] = 'r';
	//		_Filename[i+6] = 'a';
	//		_Filename[i+7] = 'w';
	//		_Filename[i+8] = '\0';
	//	}
	//}
	strcpy(_Filename,originalFilename);
}

void BinarySubimageWriter::ChangeSubimage( vtkImageData* image ){
	_image->ShallowCopy( image );
}

void BinarySubimageWriter::SetSubimageExtent( int* extent  ){
	_subimageExtent = extent;
}

void BinarySubimageWriter::SetImageDimensions ( int dimX, int dimY, int dimZ ){
	/*std::cout << "Set image dimensions" << std::endl;*/
	_imageDimensions[0] = dimX;
	_imageDimensions[1] = dimY;
	_imageDimensions[2] = dimZ;
}

void BinarySubimageWriter::WriteToFile ( ) {
	//exi : valeur initiale de l'extent de x / initial value of the x extent
	//exf : valeur finale de l'extent de y / final value of the y extent
	//Pareil pour le reste / Similar for the rest
	//A FAIRE: Les dimensions de l'image sont des fouteurs de merde. Il faut taffer dessus!
	/*std::cout << "			Image dimensions " << _imageDimensions[0] << " " << _imageDimensions[1] << " " << _imageDimensions[2] << std::endl;
	std::cout << "Pixel size " << _pixelSize << std::endl;*/

	int exi = _subimageExtent[0];
	int exf = _subimageExtent[1];

	int eyi = _subimageExtent[2];
	int eyf = _subimageExtent[3];

	int ezi = _subimageExtent[4];
	int ezf = _subimageExtent[5];
	/*std::cout << "			Extent de chaque image " << exi << " " << exf << " " << eyi << " " << eyf << " " << ezi << " " << ezf << " " << std::endl;*/

	//std::cout << "		Before writing the subvolume into the file " << std::endl;
	for( int zSV = ezi; zSV <= ezf ; zSV++ ){
		for ( int ySV = eyi; ySV <= eyf; ySV++){
			for ( int xSV = exi; xSV <= exf; xSV++ ){
				//std::cout << "X " << xSV << " Y " << ySV << " Z " << zSV << std::endl;
				char* pixel = (char*)(_image->GetScalarPointer(xSV, ySV, zSV));
				//std::cout << "Pixel value " << pixel << std::endl;
				int position = (xSV + ySV*_imageDimensions[0] + zSV*(_imageDimensions[0]*_imageDimensions[1]))*(_pixelSize)*_numberOfChannels/*_pixelSize*/;
				//std::cout << "Pixel Value " << (unsigned short*)pixel << std::endl;
 
				//_fseeki64(binary_file, position, ios::beg);
				binary_file.seekg(position, ios::beg);
				binary_file.write((char*)pixel, _pixelSize);
			}
		}
	}
	//std::cout << "		Subvolume written into the file " << std::endl;
}

//void BinarySubimageWriter::WritePixel ( ){
//	switch(_pixelType){
//	}
//}

void BinarySubimageWriter::OpenFile(){
	std::cout << "The file that is going to be opened " << _Filename << std::endl;
	binary_file.open(_Filename, fstream::out | fstream::binary );
}
void BinarySubimageWriter::CloseFile(){
	binary_file.close();
}