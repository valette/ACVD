#ifndef __vtkMHDFileWriter_h
#define __vtkMHDFileWriter_h

#include "vtkImageData.h"

class vtkMHDFileWriter {
public:
	vtkMHDFileWriter();
	~vtkMHDFileWriter();
	
	void SetMHDFileName(char* fileName);
	void SetObjectType(char* objectType);
	void SetNDimensions(int nDimensions);
	void SetDimensions(int dimX, int dimY, int dimZ);
	void SetElementType(int elementType);
	void SetHeaderSize(int size);
	void SetNumberOfChannels(int numberChn);
	void SetElementSpacing(double spacingX, double spacingY, double spacingZ);
	void SetPosition(double positionX, double positionY, double positionZ);
	void SetOffset(double offsetX, double offsetY, double offsetZ);
	void SetCenterOfRotation(double rotationX, double rotationY, double rotationZ);
	void SetIsBinaryData(bool isBinaryData);
	void SetIsByteOrderMSB(bool isByteOrderMSB);
	void SetIsCompressedData(bool isCompressedData);
	void SetElementDataFile(char* elementDataFile);	
	void WriteToFile();
private:
    char	 _mhdFileName[255];
	char	 _objectType[100];
	int		 _nDimensions;
	int		 _dimensions[3];
	char	 _elementType[200];
	int		 _headerSize;
	int		 _numberOfChannels;
	double   _elementSpacing[3];
	double   _position[3];
	double   _offset[3];
	double   _centerOfRotation[3];
	char*	 _binaryDataStr;
	char*	 _orderMSBStr;
	char*	 _compressedDataStr;
	char	 _elementDataFile[255];
	
	ofstream mhdFile;
	
};

#endif