#include "vtkMHDFileWriter.h"

vtkMHDFileWriter::vtkMHDFileWriter(){
	_nDimensions = 0;
	_headerSize = 0;
	_numberOfChannels = 0;
}

vtkMHDFileWriter::~vtkMHDFileWriter(){
}

void vtkMHDFileWriter::SetMHDFileName(char* filename){
	strcpy(_mhdFileName, filename);
}

void vtkMHDFileWriter::SetObjectType(char* objectType){
	strcpy(_objectType, objectType);
}

void vtkMHDFileWriter::SetNDimensions(int dimensions){
	_nDimensions = dimensions;
}

void vtkMHDFileWriter::SetDimensions(int dimX, int dimY, int dimZ){
	_dimensions[0] = dimX;
	_dimensions[1] = dimY;
	_dimensions[2] = dimZ;
}

void vtkMHDFileWriter::SetElementType(int type){
	switch(type)
	{
	case VTK_SIGNED_CHAR:
		strcpy(_elementType, "MET_CHAR");
		break;
	case VTK_UNSIGNED_CHAR:
		strcpy(_elementType, "MET_UCHAR");
		break;
	case VTK_SHORT:
		strcpy(_elementType, "MET_SHORT");
		break;
	case VTK_UNSIGNED_SHORT:
		strcpy(_elementType, "MET_USHORT");
		break;
	case VTK_INT:
		strcpy(_elementType, "MET_INT");
		break;
	case VTK_UNSIGNED_INT:
		strcpy(_elementType, "MET_UINT");
		break;
	case VTK_LONG:
		strcpy(_elementType, "MET_LONG");
		break;
	case VTK_UNSIGNED_LONG:
		strcpy(_elementType, "MET_ULONG");
		break;
	case VTK_FLOAT:
		strcpy(_elementType, "MET_FLOAT");
		break;
	case VTK_DOUBLE:
		strcpy(_elementType, "MET_DOUBLE");
		break;
	}
}

void vtkMHDFileWriter::SetHeaderSize(int size){
	_headerSize = size;
}

void vtkMHDFileWriter::SetNumberOfChannels(int numberOfChannels){
	_numberOfChannels = numberOfChannels;
}

void vtkMHDFileWriter::SetElementSpacing(double spacingX, double spacingY, double spacingZ){
	_elementSpacing[0] = spacingX;
	_elementSpacing[1] = spacingY;
	_elementSpacing[2] = spacingZ;
}

void vtkMHDFileWriter::SetPosition(double positionX, double positionY, double positionZ){
	_position[0] = positionX;
	_position[1] = positionY;
	_position[2] = positionZ;
}

void vtkMHDFileWriter::SetOffset(double offsetX, double offsetY, double offsetZ){
	_offset[0] = offsetX;
	_offset[1] = offsetY;
	_offset[2] = offsetZ;
}

void vtkMHDFileWriter::SetCenterOfRotation(double centerX, double centerY, double centerZ){
	_centerOfRotation[0] = centerX;
	_centerOfRotation[1] = centerY;
	_centerOfRotation[2] = centerZ;
}

void vtkMHDFileWriter::SetIsBinaryData(bool isBinaryData){
	if(isBinaryData)
		_binaryDataStr = "TRUE";
	else
		_binaryDataStr = "FALSE";
}

void vtkMHDFileWriter::SetIsByteOrderMSB(bool isByteOrderMSB){
	if(isByteOrderMSB)
		_orderMSBStr = "TRUE";
	else
		_orderMSBStr = "FALSE";
}

void vtkMHDFileWriter::SetIsCompressedData(bool isCompressedData){
	if(isCompressedData)
		_compressedDataStr = "TRUE";
	else
		_compressedDataStr = "FALSE";
}

void vtkMHDFileWriter::SetElementDataFile(char* elementDataFile){
	strcpy(_elementDataFile, elementDataFile);
}

void vtkMHDFileWriter::WriteToFile ( ){
	std::cout << " Element data file " << _elementDataFile << std::endl;
	std::cout << "	MhdFileName " << _mhdFileName << std::endl;
	std::cout << "	Begin writing the file " << std::endl;
	mhdFile.open(_mhdFileName);
	mhdFile << "ObjectType = " << _objectType << "\n";
	mhdFile << "NDims = " << _nDimensions << "\n";
	mhdFile << "DimSize = " << _dimensions[0] << " " << _dimensions[1] << " " << _dimensions[2] <<"\n";
	mhdFile << "ElementType = " << _elementType <<"\n";
	mhdFile << "HeaderSize = " << _headerSize << "\n";
	mhdFile << "ElementNumberOfChannels = " << _numberOfChannels << "\n";
	mhdFile << "ElementSpacing = " << _elementSpacing[0] << " " << _elementSpacing[1] << " " << _elementSpacing[2] << "\n";
	mhdFile << "Position = " << _position[0] << " " << _position[1] << " " << _position[2] << "\n";
	mhdFile << "Offset = " << _offset[0] << " " << _offset[1] << " " << _offset[2] << "\n";
	mhdFile << "CenterOfRotation = " << _centerOfRotation[0] << " " << _centerOfRotation[1] << " " << _centerOfRotation[2] << "\n";
	mhdFile << "BinaryData = " << _binaryDataStr <<"\n";
	mhdFile << "BinaryDataByteOrderMSB = " << _orderMSBStr << "\n";
	mhdFile << "CompressedData = " << _compressedDataStr << "\n";
	mhdFile << "ElementDataFile = " << _elementDataFile;
	mhdFile.close();
	std::cout << "	Finish writing the file " << std::endl;

}