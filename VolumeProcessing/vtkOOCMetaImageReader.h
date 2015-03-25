/*=========================================================================

  Program:     Visualization Toolkit
  Module:      $RCSfile: vtkOOCMetaImageReader.h,v $
  Authored by: Jose Vargas Casadiego, Sebastien Valette
  
  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkOOCMetaImageReader - Superclass of transformable binary file readers.
// .SECTION Description
// vtkOOCMetaImageReader provides methods needed to read a region from a mhd file.
// It supports both transforms and masks on the input data, but as a result
// is more complicated and slower than its parent class vtkImageReader2.

// .SECTION See Also
// vtkImageReader vtkMetaImageReader

#ifndef __vtkOOCMetaImageReader_h
#define __vtkOOCMetaImageReader_h
#include "vtkObjectFactory.h"
#include "vtkImageReader2.h"
//#include <windows.h>

//BTX
namespace vtkmetaio { class MetaImage; } // forward declaration  
//ETX

class vtkTransform;

#define VTK_FILE_BYTE_ORDER_BIG_ENDIAN 0
#define VTK_FILE_BYTE_ORDER_LITTLE_ENDIAN 1

class VTK_EXPORT vtkOOCMetaImageReader : public vtkImageReader2
{
public:
	static vtkOOCMetaImageReader *New();
	//{return new vtkOOCMetaImageReader;};

	vtkGetMacro(XMin, int);
	vtkSetMacro(XMin, int);
	vtkGetMacro(XMax, int);
	vtkSetMacro(XMax, int);
	vtkGetMacro(YMin, int);
	vtkSetMacro(YMin, int);
	vtkGetMacro(YMax, int);
	vtkSetMacro(YMax, int);
	vtkGetMacro(ZMin, int);
	vtkSetMacro(ZMin, int);
	vtkGetMacro(ZMax, int);
	vtkSetMacro(ZMax, int);



  //Meta-File reading functionality
    virtual const char * GetFileExtensions() 
    { return ".mhd .mha"; }

  virtual const char * GetDescriptiveName() 
    { return "MetaIO Library: MetaImage"; }
  //Meta-file reading functionality 

  /**
  ** Additional methods required for extracting all the information that's contained within
  ** .mhd files
  **/
  void ReadImageParameters();
  double * GetPixelSpacing()
    { return this->GetDataSpacing(); }
  int GetWidth()
    { return (this->GetDataExtent()[1] - this->GetDataExtent()[0] + 1); }
  int GetHeight()
    { return (this->GetDataExtent()[3] - this->GetDataExtent()[2] + 1); }
  double * GetImagePositionPatient()
    { return this->GetDataOrigin(); }
  int GetNumberOfComponents()
    { return this->GetNumberOfScalarComponents(); }
  int GetPixelRepresentation()
    { return this->GetDataScalarType(); }
  int GetDataByteOrder(); 

  /**
  * Additional macros required for extracting all the information that's contained within .mhd files
  **/
  vtkGetMacro(RescaleSlope, double);
  vtkGetMacro(RescaleOffset, double);
  vtkGetMacro(BitsAllocated, int);
  vtkGetStringMacro(DistanceUnits);
  vtkGetStringMacro(AnatomicalOrientation);
  vtkGetMacro(GantryAngle, double); 
  vtkGetStringMacro(PatientName);
  vtkGetStringMacro(PatientID);
  vtkGetStringMacro(Date);
  vtkGetStringMacro(Series);
  vtkGetStringMacro(ImageNumber);
  vtkGetStringMacro(Modality);
  vtkGetStringMacro(StudyID);
  vtkGetStringMacro(StudyUID);
  vtkGetStringMacro(TransferSyntaxUID);

  // Description:
  // Test whether the file with the given name can be read by this
  // reader.
  virtual int CanReadFile(const char* name);
  
  //  vtkTypeRevisionMacro(vtkOOCMetaImageReader,vtkImageReader2);
  void PrintSelf(ostream& os, vtkIndent indent);   

  // Description:
  // Set/get the data VOI. You can limit the reader to only
  // read a subset of the data. 
  vtkSetVector6Macro(DataVOI,int);
  vtkGetVector6Macro(DataVOI,int);
  
  // Description:
  // Set/Get the Data mask.
  vtkGetMacro(DataMask,unsigned short);
  void SetDataMask(int val) 
    {
      if (val == this->DataMask)
        {
        return;
        }
      this->DataMask = static_cast<unsigned short>(val);
      this->Modified();
    }
  
  // Description:
  // Set/Get transformation matrix to transform the data from slice space
  // into world space. This matrix must be a permutation matrix. To qualify,
  // the sums of the rows must be + or - 1.
  virtual void SetTransform(vtkTransform*);
  vtkGetObjectMacro(Transform,vtkTransform);

  // Warning !!!
  // following should only be used by methods or template helpers, not users
  void ComputeInverseTransformedExtent(int inExtent[6],
                                       int outExtent[6]);
  void ComputeInverseTransformedIncrements(vtkIdType inIncr[3],
                                           vtkIdType outIncr[3]);

  int OpenAndSeekFile(int extent[6], int slice);
  
  // Description:
  // Set/get the scalar array name for this data set.
  vtkSetStringMacro(ScalarArrayName);
  vtkGetStringMacro(ScalarArrayName);
  
protected:
  vtkOOCMetaImageReader();
  ~vtkOOCMetaImageReader();

  unsigned short DataMask;  // Mask each pixel with ...

  vtkTransform *Transform;

  void ComputeTransformedSpacing (double Spacing[3]);
  void ComputeTransformedOrigin (double origin[3]);
  void ComputeTransformedExtent(int inExtent[6],
                                int outExtent[6]);
  void ComputeTransformedIncrements(vtkIdType inIncr[3],
                                    vtkIdType outIncr[3]);

  int DataVOI[6];
  int XMin, XMax, YMin, YMax, ZMin, ZMax;
  
  char *ScalarArrayName;
  
  //This method is not overloaded in vtkImageReader; it is, however, in vtkMetaImageReader (quite predictably,
  //it must be said). Hence, it shall replicate the latter method. 
  void ExecuteInformation();
  
  //Here, the .mhd file information is requested by the reader.
  //As such, this method has to duplicate the behavior of the vtkMetaImageReader
  //Nonetheless, it is important to note that the implementations contained within
  //vtkImageReader and vtkMetaImageReader difer by more than small details. As such,
  //a middle-ground must be found between both
  virtual int RequestInformation(vtkInformation* request,
                                 vtkInformationVector** inputVector,
                                 vtkInformationVector* outputVector);

	
  //This is the method where the actual reading of the file takes places. As such, it should
  //replicate vtkImageReader's functionality because it's this class that permets the lecture of
  //regions of interest
  void ExecuteData(vtkDataObject *data, vtkInformation* outInfo);
private:
  vtkOOCMetaImageReader(const vtkOOCMetaImageReader&);  // Not implemented.
  void operator=(const vtkOOCMetaImageReader&);  // Not implemented.

  //BTX
  vtkmetaio::MetaImage *MetaImagePtr;
  //ETX

  double GantryAngle;
  char PatientName[255];
  char PatientID[255];
  char Date[255];
  char Series[255];
  char Study[255];
  char ImageNumber[255];
  char Modality[255];
  char StudyID[255];
  char StudyUID[255];
  char TransferSyntaxUID[255];

  double RescaleSlope;
  double RescaleOffset;
  int BitsAllocated;
  char DistanceUnits[255];
  char AnatomicalOrientation[255];

  bool parametersAlreadyRead;
};

#endif

