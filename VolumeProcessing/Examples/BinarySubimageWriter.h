#ifndef __BinarySubimageWriter_h
#define __BinarySubimageWriter_h

#include "vtkImageData.h"

class BinarySubimageWriter {
	private:

		fstream                       binary_file;
		static BinarySubimageWriter   *instance;
		int                           _numberOfChannels;
		int						      _pixelSize;
		int                           _pixelType;
		char                          _Filename[300];
		int*                          _subimageExtent;
		int*                          _imageDimensions;
		vtkImageData*                 _image;
	public:
		BinarySubimageWriter();
		~BinarySubimageWriter();
		/*static BinarySubimageWriter* GetInstance();*/
		void   SetNumberOfChannels ( int numberOfChannels );
		void   SetPixelSizeAndType ( int pixelType );
		void   AssignFilename ( char* originalFilename );
		void   SetSubimageExtent( int* extent );
		void   SetImageDimensions ( int dimX, int dimY, int dimZ );
		void   ChangeSubimage( vtkImageData* image );
		void   WriteToFile();
		void   CloseFile();
		void   OpenFile();
};
//BinarySubimageWriter* BinarySubimageWriter::instance = NULL;

#endif