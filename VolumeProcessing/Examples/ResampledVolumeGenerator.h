
#ifndef __ResampledVolumeGenerator_h
#define __ResampledVolumeGenerator_h

#include "vtkOOCMetaImageReader.h"
#include "vtkImageResample.h"
#include "vtkImageData.h"
#include "BinarySubimageWriter.h"
#include "vtkMHDFileWriter.h"
#include "vtkImageAppendComponents.h"
#include "vtkImageAppend.h"


class ResampledVolumeGenerator{
	private:
		bool _write;

		int _divisionFactorX;
		int _divisionFactorY;
		int _divisionFactorZ;
		
		double _resamplingFactorX;
		double _resamplingFactorY;
		double _resamplingFactorZ;
		
		int _dimensionX;
		int _dimensionY;
		int _dimensionZ;

		int _resampledDimensionX;
		int _resampledDimensionY;
		int _resampledDimensionZ;
		
		int _stepX;
		int _stepY;
		int _stepZ;
		
		int _dataType;
		int _dataTypeSize;
		
		char _mhdFilename[300];
		char _rawFilename[300];
		char _filename[300];
		char _originalFilename[300];

		char* _truncatedFilename;
		char  _elementDataName[300];

		void  ConfigureFilenames ( );
		void GetImageInformation();
		void ObtainTruncatedFilename ( );
		void RawDataFilename();
		void MHDFilename();
		void AssignElementDataFile();

		BinarySubimageWriter    *writer;
		vtkOOCMetaImageReader   *metaReader;
		vtkImageResample        *resampler;
		vtkMHDFileWriter        *mhdwriter;
		vtkImageData            *inter;
		
	protected:
		virtual void StartWork();

	public:
		ResampledVolumeGenerator ( );
		~ResampledVolumeGenerator ( );
		
		char* GetNewMHDFileName ( );
		void  SetWriteToRawFile ( bool write );
		void  SetDivisionFactors ( int divisionFactorX, int divisionFactorY, int divisionFactorZ );
		void  SetResamplingFactors ( double resamplingFactorX, double resamplingFactorY, double resamplingFactorZ );
		void  SetFileName ( char* filename );
		void  UnthreadedExecute ( );
	
};
#endif
