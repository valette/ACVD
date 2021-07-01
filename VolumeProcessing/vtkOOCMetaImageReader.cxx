/*=========================================================================

Program:   Visualization Toolkit
Module:    $RCSfile: vtkOOCMetaImageReader.cxx,v $

Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
All rights reserved.
See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkOOCMetaImageReader.h"

#include <cmath>

#include <vtkByteSwap.h>
#include <vtkDataArray.h>
#include <vtkImageData.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>

#include <vtkPointData.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkTransform.h>

#include <vtkmetaio/metaTypes.h>
#include <vtkmetaio/metaUtils.h>
#include <vtkmetaio/metaEvent.h>
#include <vtkmetaio/metaObject.h>
#include <vtkmetaio/metaImageTypes.h>
#include <vtkmetaio/metaImageUtils.h>
#include <vtkmetaio/metaImage.h>

vtkStandardNewMacro(vtkOOCMetaImageReader);

//vtkCxxRevisionMacro(vtkOOCMetaImageReader, "$Revision: 1.1 $");

vtkCxxSetObjectMacro(vtkOOCMetaImageReader,Transform,vtkTransform);

#ifdef read
#undef read
#endif

#ifdef close
#undef close
#endif

//----------------------------------------------------------------------------
vtkOOCMetaImageReader::vtkOOCMetaImageReader()
{

	//MetaImageReader
	GantryAngle = 0;
	strcpy(PatientName, "?");
	strcpy(PatientID, "?");
	strcpy(Date, "?");
	strcpy(Series, "?");
	strcpy(Study, "?");
	strcpy(ImageNumber, "?");
	strcpy(Modality, "?");
	strcpy(StudyID, "?");
	strcpy(StudyUID, "?");
	strcpy(TransferSyntaxUID, "?");

	RescaleSlope = 1;
	RescaleOffset = 0;
	BitsAllocated = 0;
	strcpy(DistanceUnits, "mm");
	strcpy(AnatomicalOrientation, "RAS");
	this->MetaImagePtr = new vtkmetaio::MetaImage;
	this->FileLowerLeft = 1;

	//ImageReader
	int idx;

	for (int i = 0; i < 6; ++i) {
		this->DataVOI[i] = VTK_INT_MIN;
	}

	this->XMin = this->XMax = this->YMin = this->YMax = this->ZMin = this->ZMax = VTK_INT_MIN;

	// Left over from short reader
	this->DataMask = 0xffff;
	this->Transform = NULL;

	this->ScalarArrayName = NULL;
	this->SetScalarArrayName("MetaImage");

	parametersAlreadyRead = false;
}

//----------------------------------------------------------------------------
vtkOOCMetaImageReader::~vtkOOCMetaImageReader()
{ 
	delete this->MetaImagePtr;
	this->SetTransform(NULL);
	this->SetScalarArrayName(NULL);
}

//----------------------------------------------------------------------------
void vtkOOCMetaImageReader::PrintSelf(ostream& os, vtkIndent indent)
{
	int idx;

	this->Superclass::PrintSelf(os,indent);

	os << indent << "Data Mask: " << this->DataMask << "\n";
	os << indent << "DataVOI: (" << this->DataVOI[0];
	for (idx = 1; idx < 6; ++idx)
	{
		os << ", " << this->DataVOI[idx];
	}
	os << ")\n";
	if ( this->Transform )
	{
		os << indent << "Transform: " << this->Transform << "\n";
	}
	else
	{
		os << indent << "Transform: (none)\n";
	}

	os << indent << "ScalarArrayName: "
		<< (this->ScalarArrayName ? this->ScalarArrayName : "(none)") << endl;
}

void vtkOOCMetaImageReader::ExecuteInformation()
{
	if(!parametersAlreadyRead){
		ReadImageParameters ( );
	}
}

void vtkOOCMetaImageReader::ReadImageParameters ( ) {
	if(!this->FileName)
	{
		vtkErrorMacro( << "A filename was not specified." );
		return;
	}

	if(!this->MetaImagePtr->Read(this->FileName, false))
	{
		vtkErrorMacro( << "MetaImage cannot parse file." );
		return;
	}

	if (this->MetaImagePtr->BinaryDataByteOrderMSB()) {
		this->SetDataByteOrderToBigEndian();
	}

	this->SetFileDimensionality(this->MetaImagePtr->NDims());
	if ( FileDimensionality <= 0 || FileDimensionality >= 4)
	{
		vtkErrorMacro(
			<< "Only understands image data of 1, 2, and 3 dimensions. "
			<< "This image has " << FileDimensionality << " dimensions");
		return;
	}

	//The cout qui montre la dimensionalité du fichier
	vtkDebugMacro(<< "* This image has " << FileDimensionality << " dimensions");

	int i;

	switch(this->MetaImagePtr->ElementType())
	{
	default:
		vtkErrorMacro(<< "Unknown data type: " 
			<< this->MetaImagePtr->ElementType());
		return;
	case vtkmetaio::MET_CHAR:
	case vtkmetaio::MET_CHAR_ARRAY:
		this->DataScalarType = VTK_SIGNED_CHAR;
		break;
	case vtkmetaio::MET_UCHAR:
	case vtkmetaio::MET_UCHAR_ARRAY:
		this->DataScalarType = VTK_UNSIGNED_CHAR;
		break;
	case vtkmetaio::MET_SHORT:
	case vtkmetaio::MET_SHORT_ARRAY:
		this->DataScalarType = VTK_SHORT;
		break;
	case vtkmetaio::MET_USHORT:
	case vtkmetaio::MET_USHORT_ARRAY:
		this->DataScalarType = VTK_UNSIGNED_SHORT;
		break;
	case vtkmetaio::MET_INT:
	case vtkmetaio::MET_INT_ARRAY:
		this->DataScalarType = VTK_INT;
		break;
	case vtkmetaio::MET_UINT:
	case vtkmetaio::MET_UINT_ARRAY:
		this->DataScalarType = VTK_UNSIGNED_INT;
		break;
	case vtkmetaio::MET_LONG:	
	case vtkmetaio::MET_LONG_ARRAY:
		this->DataScalarType = VTK_LONG;
		break;
	case vtkmetaio::MET_ULONG:
	case vtkmetaio::MET_ULONG_ARRAY:
		this->DataScalarType = VTK_UNSIGNED_LONG;
		break;
	case vtkmetaio::MET_FLOAT:
		this->DataScalarType = VTK_FLOAT;
		break;
	case vtkmetaio::MET_DOUBLE:
		this->DataScalarType = VTK_DOUBLE;
		break;
	}

	int extent[6]={0,0,0,0,0,0};
	double spacing[3]={1.0, 1.0, 1.0};
	double origin[3]={0.0, 0.0, 0.0};
	for(i=0; i<FileDimensionality; i++)
	{
		extent[2*i] = 0;
		extent[2*i+1] = this->MetaImagePtr->DimSize(i)-1;
		spacing[i] = fabs(this->MetaImagePtr->ElementSpacing(i));
		origin[i] = this->MetaImagePtr->Position(i);
	}

	this->SetNumberOfScalarComponents(
		this->MetaImagePtr->ElementNumberOfChannels());
	this->SetDataExtent(extent);
	this->SetDataSpacing(spacing);
	this->SetDataOrigin(origin);
	this->SetHeaderSize(this->MetaImagePtr->HeaderSize());
	this->FileLowerLeftOn();

	//Do you have to specify the ordering of the bytes ?
	/*if(this->MetaImagePtr->BinaryDataByteOrderMSB())
	this->SetDataByteOrder(1);
	else
	this->SetDataByteOrder(0);*/

	switch(this->MetaImagePtr->DistanceUnits())
	{
	default:
	case vtkmetaio::MET_DISTANCE_UNITS_UM:
		{
			strcpy(DistanceUnits, "um");
			break;
		}
	case vtkmetaio::MET_DISTANCE_UNITS_MM:
		{
			strcpy(DistanceUnits, "mm");
			break;
		}
	case vtkmetaio::MET_DISTANCE_UNITS_CM:
		{
			strcpy(DistanceUnits, "cm");
			break;
		}
	}

	strcpy(AnatomicalOrientation,
		this->MetaImagePtr->AnatomicalOrientationAcronym());

	vtkmetaio::MET_SizeOfType(this->MetaImagePtr->ElementType(), &BitsAllocated);

	RescaleSlope = this->MetaImagePtr->ElementToIntensityFunctionSlope();
	RescaleOffset = this->MetaImagePtr->ElementToIntensityFunctionOffset();

	if(this->MetaImagePtr->Modality() == vtkmetaio::MET_MOD_CT)
	{
		strcpy(Modality, "CT");
	}
	else if(this->MetaImagePtr->Modality() == vtkmetaio::MET_MOD_MR)
	{
		strcpy(Modality, "MR");
	}
	else
	{
		strcpy(Modality, "?");
	}

	bool usePath;
	char pathName[255];
	char fName[255];

	cout << "FIXME!!! vtkOOCMetaImageREader.cxx" << endl;
	exit( 1 );/*
	usePath = vtkmetaio::MET_GetFilePath( FileName, pathName );

	if(usePath)
	{
		sprintf(fName, "%s%s", pathName, MetaImagePtr->ElementDataFileName());
	}
	else
	{
		strcpy(fName, MetaImagePtr->ElementDataFileName());
	}
	strcpy( FileName, fName );
	parametersAlreadyRead = true;*/
}

// This method returns the largest data that can be generated.
int vtkOOCMetaImageReader::RequestInformation (
	vtkInformation       * vtkNotUsed( request ),
	vtkInformationVector** vtkNotUsed( inputVector ),
	vtkInformationVector * outputVector)
{

	// call the old method to help with backwards compatiblity
	this->ExecuteInformation();

	// get the info objects
	vtkInformation* outInfo = outputVector->GetInformationObject(0);
	double spacing[3];
	int extent[6];
	double origin[3];

	// set the extent, if the VOI has not been set then default to
	// the DataExtent
	int VOI[6];
	for (int i = 0; i < 6; i++) {
		if (this->DataVOI[i] != VTK_INT_MIN) {
			VOI[i] = this->DataVOI[i];
		} else {
			VOI[i] = this->DataExtent[i];
		}
	}

	if (this->XMin != VTK_INT_MIN) {
		VOI[0] = this->XMin;
	}
	if (this->XMax != VTK_INT_MIN) {
		VOI[1] = this->XMax;
	}
	if (this->YMin != VTK_INT_MIN) {
		VOI[2] = this->YMin;
	}
	if (this->YMax != VTK_INT_MIN) {
		VOI[3] = this->YMax;
	}
	if (this->ZMin != VTK_INT_MIN) {
		VOI[4] = this->ZMin;
	}
	if (this->ZMax != VTK_INT_MIN) {
		VOI[5] = this->ZMax;
	}

	this->ComputeTransformedExtent(VOI, extent);

	outInfo->Set(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(),extent, 6);

	this->ComputeTransformedSpacing(spacing);
	outInfo->Set(vtkDataObject::SPACING(), this->DataSpacing, 3);

	this->ComputeTransformedOrigin(origin);
	outInfo->Set(vtkDataObject::ORIGIN(),  this->DataOrigin, 3);

	vtkDataObject::SetPointDataActiveScalarInfo(outInfo, this->DataScalarType,
		this->NumberOfScalarComponents);
	return 1;


}

int vtkOOCMetaImageReader::OpenAndSeekFile(int dataExtent[6], int idx)
{
	long long streamStart;

	if (!this->FileName && !this->FilePattern)
	{
		vtkErrorMacro(<<"Either a FileName or FilePattern must be specified.");
		return 0;
	}
	this->ComputeInternalFileName(idx);
	this->OpenFile();
	if ( !this->File )
	{
		return 0;
	}
	// convert data extent into constants that can be used to seek.
	streamStart = 
		(dataExtent[0] - this->DataExtent[0]) * this->DataIncrements[0];

	if (this->FileLowerLeft)
	{
		streamStart = streamStart + 
			(dataExtent[2] - this->DataExtent[2]) * this->DataIncrements[1];
	}
	else
	{
		streamStart = streamStart + 
			(this->DataExtent[3] - this->DataExtent[2] - dataExtent[2]) * 
			this->DataIncrements[1];
	}

	// handle three and four dimensional files
	if (this->GetFileDimensionality() >= 3)
	{
		streamStart = streamStart + 
			(dataExtent[4] - this->DataExtent[4]) * this->DataIncrements[2];
	}

	streamStart += this->GetHeaderSize(idx);

	// error checking
	//JFVC modif
	this->File->seekg(static_cast<long long>(streamStart), ios::beg);
	if (this->File->fail())
	{
		vtkErrorMacro(<< "File operation failed: " << streamStart << ", ext: "
			<< dataExtent[0] << ", " << dataExtent[1] << ", "
			<< dataExtent[2] << ", " << dataExtent[3] << ", "
			<< dataExtent[4] << ", " << dataExtent[5]);
		vtkErrorMacro(<< "Header size: " << this->GetHeaderSize(idx) << ", file ext: "
			<< this->DataExtent[0] << ", " << this->DataExtent[1] << ", "
			<< this->DataExtent[2] << ", " << this->DataExtent[3] << ", "
			<< this->DataExtent[4] << ", " << this->DataExtent[5]);
		return 0;
	}
	return 1;        
}

//----------------------------------------------------------------------------
// This function reads in one data of data.
// templated to handle different data types.
template <class IT, class OT>
void vtkOOCMetaImageReaderUpdate2(vtkOOCMetaImageReader *self, vtkImageData *data,
								  IT *inPtr, OT *outPtr)
{
	vtkIdType inIncr[3], outIncr[3];
	OT *outPtr0, *outPtr1, *outPtr2;
	long long streamSkip0, streamSkip1;
	unsigned long streamRead;
	int idx0, idx1, idx2, pixelRead;
	unsigned char *buf;
	int inExtent[6];
	int dataExtent[6];
	int comp, pixelSkip;
	long long filePos, correction = 0;
	unsigned long count = 0;
	unsigned short DataMask;
	unsigned long target;

	// Get the requested extents.
	data->GetExtent(inExtent);
	// Convert them into to the extent needed from the file. 
	self->ComputeInverseTransformedExtent(inExtent,dataExtent);

	// get and transform the increments
	data->GetIncrements(inIncr);
	self->ComputeInverseTransformedIncrements(inIncr,outIncr);

	DataMask = self->GetDataMask();

	// compute outPtr2 
	outPtr2 = outPtr;
	if (outIncr[0] < 0) 
	{
		outPtr2 = outPtr2 - outIncr[0]*(dataExtent[1] - dataExtent[0]);
	}
	if (outIncr[1] < 0) 
	{
		outPtr2 = outPtr2 - outIncr[1]*(dataExtent[3] - dataExtent[2]);
	}
	if (outIncr[2] < 0) 
	{
		outPtr2 = outPtr2 - outIncr[2]*(dataExtent[5] - dataExtent[4]);
	}

	// length of a row, num pixels read at a time
	pixelRead = dataExtent[1] - dataExtent[0] + 1; 
	streamRead = static_cast<unsigned long>(pixelRead * 
		self->GetDataIncrements()[0]);  
	streamSkip0 = (long)(self->GetDataIncrements()[1] - streamRead);
	streamSkip1 = (long)(self->GetDataIncrements()[2] - 
		(dataExtent[3] - dataExtent[2] + 1)* self->GetDataIncrements()[1]);
	pixelSkip = data->GetNumberOfScalarComponents();

	// read from the bottom up
	if (!self->GetFileLowerLeft()) 
	{
		streamSkip0 = (long)(-static_cast<long>(streamRead) 
			- self->GetDataIncrements()[1]);
		streamSkip1 = (long)(self->GetDataIncrements()[2] + 
			(dataExtent[3] - dataExtent[2] + 1)* self->GetDataIncrements()[1]);
	}


	// create a buffer to hold a row of the data
	buf = new unsigned char[streamRead];

	target = (unsigned long)((dataExtent[5]-dataExtent[4]+1)*
		(dataExtent[3]-dataExtent[2]+1)/50.0);
	target++;

	// read the data row by row
	if (self->GetFileDimensionality() == 3)
	{
		if ( !self->OpenAndSeekFile(dataExtent,0) )
		{
			delete [] buf;
			return;
		}
	}
	for (idx2 = dataExtent[4]; idx2 <= dataExtent[5]; ++idx2)
	{
		if (self->GetFileDimensionality() == 2)
		{
			if ( !self->OpenAndSeekFile(dataExtent,idx2) )
			{
				delete [] buf;
				return;
			}
		}
		outPtr1 = outPtr2;
		for (idx1 = dataExtent[2]; 
			!self->AbortExecute && idx1 <= dataExtent[3]; ++idx1)
		{
			if (!(count%target))
			{
				self->UpdateProgress(count/(50.0*target));
			}
			count++;
			outPtr0 = outPtr1;

			// read the row.
			self->GetFile()->read((char *)buf, streamRead);
#ifdef __APPLE_CC__
			if (static_cast<unsigned long>(self->GetFile()->gcount()) != streamRead)
				// Apple's gcc3 returns fail when reading _to_ eof
#else
			if ( static_cast<unsigned long>(self->GetFile()->gcount()) != 
				streamRead || self->GetFile()->fail())
#endif
			{
				vtkGenericWarningMacro("File operation failed. row = " << idx1
					<< ", Tried to Read = " << streamRead
					<< ", Read = " << self->GetFile()->gcount()
					<< ", Skip0 = " << streamSkip0
					<< ", Skip1 = " << streamSkip1
					<< ", FilePos = " << static_cast<long long>(self->GetFile()->tellg()));
				delete [] buf;
				return;
			}
			// handle swapping
			if (self->GetSwapBytes())
			{
				// pixelSkip is the number of components in data
				vtkByteSwap::SwapVoidRange(buf, pixelRead*pixelSkip, sizeof(IT));
			}

			// copy the bytes into the typed data
			inPtr = (IT *)(buf);
			for (idx0 = dataExtent[0]; idx0 <= dataExtent[1]; ++idx0)
			{
				// Copy pixel into the output.
				if (DataMask == 0xffff)
				{
					for (comp = 0; comp < pixelSkip; comp++)
					{
						outPtr0[comp] = (OT)(inPtr[comp]);
					}
				}
				else
				{
					// left over from short reader (what about other types.
					for (comp = 0; comp < pixelSkip; comp++)
					{
						outPtr0[comp] = (OT)((short)(inPtr[comp]) & DataMask);
					}
				}
				// move to next pixel
				inPtr += pixelSkip;
				outPtr0 += outIncr[0];
			}

			// move to the next row in the file and data
			filePos = self->GetFile()->tellg();

			/* Unfortunately this doesn't work as a fix, but I'll leave it here for a bit
			to provoke ideas.
			if (filePos == -1)
			{
			self->GetFile()->clear(self->GetFile()->rdstate() & ~ios::eofbit);
			self->GetFile()->clear(self->GetFile()->rdstate() & ~ios::failbit);
			filePos = self->GetFile()->tellg();
			}
			*/
#if defined(VTK_USE_ANSI_STDLIB ) && ((defined(__sgi) && !defined(__GNUC__)) \
	|| (__BORLANDC__>=0x0560) || defined (__APPLE_CC__))
			// this check is required for SGI's when vtk is build with VTK_USE_ANSI_STDLIB
			// seems that after a read that just reaches EOF, tellg reports a -1.
			// clear() does not work, so we have to reopen the file.
			// NB: Also Borland CBuilder 6 suffers from same trouble
			// As does Apple's gcc3
			if (filePos == -1)
			{
				self->OpenFile();
				self->GetFile()->seekg(0,ios::end);
				filePos = self->GetFile()->tellg();
			}
#endif
			// watch for case where we might rewind too much
			// if that happens, store the value in correction and apply later
			if (filePos + streamSkip0 >= 0)
			{
				self->GetFile()->seekg(static_cast<long long>(self->GetFile()->tellg()) + streamSkip0, ios::beg);
				correction = 0;
			}
			else
			{
				correction = streamSkip0;
			}
			outPtr1 += outIncr[1];
		}
		// move to the next image in the file and data
		//JF//VC : Modif
			//SetFilePointerEx( self->GetFile(),static_cast<_int64>(self->GetFile()->tellg()) + streamSkip1 + correction, 0, FILE_BEGIN);
		self->GetFile()->seekg(static_cast<long long>(self->GetFile()->tellg()) + streamSkip1 + correction, 
			ios::beg);
		outPtr2 += outIncr[2];
	}

	// delete the temporary buffer
	delete [] buf;
}


//----------------------------------------------------------------------------
// This function reads in one data of one slice.
// templated to handle different data types.
template <class T>
void vtkOOCMetaImageReaderUpdate1(vtkOOCMetaImageReader *self, vtkImageData *data, T *inPtr)
{
	void *outPtr;

	// Call the correct templated function for the input
	outPtr = data->GetScalarPointer();
	switch (data->GetScalarType())
	{
		vtkTemplateMacro(vtkOOCMetaImageReaderUpdate2(self, data, inPtr, 
			(VTK_TT *)(outPtr)));
	default:
		vtkGenericWarningMacro("Update1: Unknown data type\n");
	}  
}
//----------------------------------------------------------------------------
// This function reads a data from a file.  The datas extent/axes
// are assumed to be the same as the file extent/order.
void vtkOOCMetaImageReader::ExecuteData(vtkDataObject *output,
                                     vtkInformation* outInfo)
{
	vtkImageData *data = this->AllocateOutputData(output, outInfo);

	void *ptr = NULL;
	int *ext;

	if (!this->FileName && !this->FilePattern)
	{
		vtkErrorMacro("Either a valid FileName or FilePattern must be specified.");
		return;
	}

	ext = data->GetExtent();
	if (!data->GetPointData()->GetScalars())
	{
		return;
	}
	data->GetPointData()->GetScalars()->SetName("MetaImage");

	vtkDebugMacro("	ExecuteData, Reading extent: " << ext[0] << ", " << ext[1] << ", " 
		<< ext[2] << ", " << ext[3] << ", " << ext[4] << ", " << ext[5]);

	this->ComputeDataIncrements();

	// Call the correct templated function for the output
	switch (this->GetDataScalarType())
	{
		vtkTemplateMacro(vtkOOCMetaImageReaderUpdate1(this, data, (VTK_TT *)(ptr)));
	default:
		vtkErrorMacro(<< "UpdateFromFile: Unknown data type");
	}   
}





void vtkOOCMetaImageReader::ComputeTransformedSpacing (double Spacing[3])
{
	if (!this->Transform)
	{
		memcpy (Spacing, this->DataSpacing, 3 * sizeof (double));
	}
	else
	{
		double transformedSpacing[3];
		memcpy (transformedSpacing, this->DataSpacing, 3 * sizeof (double));
		this->Transform->TransformVector(transformedSpacing, transformedSpacing);

		for (int i = 0; i < 3; i++)
		{
			Spacing[i] = fabs(transformedSpacing[i]);
		}
		vtkDebugMacro("Transformed Spacing " << Spacing[0] << ", " << Spacing[1] << ", " << Spacing[2]);
	}
}

// if the spacing is negative we need to tranlate the origin
// basically O' = O + spacing*(dim-1) for any axis that would
// have a negative spaing
void vtkOOCMetaImageReader::ComputeTransformedOrigin (double origin[3])
{
	if (!this->Transform)
	{
		memcpy (origin, this->DataOrigin, 3 * sizeof (double));
	}
	else
	{
		double transformedOrigin[3];
		double transformedSpacing[3];
		int transformedExtent[6];

		memcpy (transformedSpacing, this->DataSpacing, 3 * sizeof (double));
		this->Transform->TransformVector(transformedSpacing, transformedSpacing);

		memcpy (transformedOrigin, this->DataOrigin, 3 * sizeof (double));
		this->Transform->TransformPoint(transformedOrigin, transformedOrigin);

		this->ComputeTransformedExtent(this->DataExtent,transformedExtent);

		for (int i = 0; i < 3; i++) 
		{
			if (transformedSpacing[i] < 0)
			{
				origin[i] = transformedOrigin[i] + transformedSpacing[i]*
					(transformedExtent[i*2+1] -  transformedExtent[i*2] + 1);
			}
			else
			{
				origin[i] = transformedOrigin[i];
			}
		}
		vtkDebugMacro("	ComputeTransformedOrigin, Transformed Origin " << origin[0] << ", " << origin[1] << ", " << origin[2]);
	}
}

void vtkOOCMetaImageReader::ComputeTransformedExtent(int inExtent[6],
													 int outExtent[6])
{
	double transformedExtent[3];
	int temp;
	int idx;
	int dataExtent[6];

	if (!this->Transform)
	{
		memcpy (outExtent, inExtent, 6 * sizeof (int));
		memcpy (dataExtent, this->DataExtent, 6 * sizeof(int));
	}
	else
	{
		// need to know how far to translate to start at 000
		// first transform the data extent
		transformedExtent[0] = this->DataExtent[0];
		transformedExtent[1] = this->DataExtent[2];
		transformedExtent[2] = this->DataExtent[4];
		this->Transform->TransformPoint(transformedExtent, transformedExtent);
		dataExtent[0] = (int) transformedExtent[0];
		dataExtent[2] = (int) transformedExtent[1];
		dataExtent[4] = (int) transformedExtent[2];

		transformedExtent[0] = this->DataExtent[1];
		transformedExtent[1] = this->DataExtent[3];
		transformedExtent[2] = this->DataExtent[5];
		this->Transform->TransformPoint(transformedExtent, transformedExtent);
		dataExtent[1] = (int) transformedExtent[0];
		dataExtent[3] = (int) transformedExtent[1];
		dataExtent[5] = (int) transformedExtent[2];

		for (idx = 0; idx < 6; idx += 2)
		{
			if (dataExtent[idx] > dataExtent[idx+1]) 
			{
				temp = dataExtent[idx];
				dataExtent[idx] = dataExtent[idx+1];
				dataExtent[idx+1] = temp;
			}
		}

		// now transform the inExtent
		transformedExtent[0] = inExtent[0];
		transformedExtent[1] = inExtent[2];
		transformedExtent[2] = inExtent[4];
		this->Transform->TransformPoint(transformedExtent, transformedExtent);
		outExtent[0] = (int) transformedExtent[0];
		outExtent[2] = (int) transformedExtent[1];
		outExtent[4] = (int) transformedExtent[2];

		transformedExtent[0] = inExtent[1];
		transformedExtent[1] = inExtent[3];
		transformedExtent[2] = inExtent[5];
		this->Transform->TransformPoint(transformedExtent, transformedExtent);
		outExtent[1] = (int) transformedExtent[0];
		outExtent[3] = (int) transformedExtent[1];
		outExtent[5] = (int) transformedExtent[2];
	}

	for (idx = 0; idx < 6; idx += 2)
	{
		if (outExtent[idx] > outExtent[idx+1]) 
		{
			temp = outExtent[idx];
			outExtent[idx] = outExtent[idx+1];
			outExtent[idx+1] = temp;
		}
		// do the slide to 000 origin by subtracting the minimum extent
		outExtent[idx] -= dataExtent[idx];
		outExtent[idx+1] -= dataExtent[idx];
	}

	vtkDebugMacro(<< "Transformed extent are:" 
		<< outExtent[0] << ", " << outExtent[1] << ", "
		<< outExtent[2] << ", " << outExtent[3] << ", "
		<< outExtent[4] << ", " << outExtent[5]);
}

void vtkOOCMetaImageReader::ComputeInverseTransformedExtent(int inExtent[6],
															int outExtent[6])
{
	double transformedExtent[3];
	int temp;
	int idx;

	if (!this->Transform)
	{
		memcpy (outExtent, inExtent, 6 * sizeof (int));
		for (idx = 0; idx < 6; idx += 2)
		{
			// do the slide to 000 origin by subtracting the minimum extent
			outExtent[idx] += this->DataExtent[idx];
			outExtent[idx+1] += this->DataExtent[idx];
		}
	}
	else
	{
		// need to know how far to translate to start at 000
		int dataExtent[6];
		// first transform the data extent
		transformedExtent[0] = this->DataExtent[0];
		transformedExtent[1] = this->DataExtent[2];
		transformedExtent[2] = this->DataExtent[4];
		this->Transform->TransformPoint(transformedExtent, transformedExtent);
		dataExtent[0] = (int) transformedExtent[0];
		dataExtent[2] = (int) transformedExtent[1];
		dataExtent[4] = (int) transformedExtent[2];

		transformedExtent[0] = this->DataExtent[1];
		transformedExtent[1] = this->DataExtent[3];
		transformedExtent[2] = this->DataExtent[5];
		this->Transform->TransformPoint(transformedExtent, transformedExtent);
		dataExtent[1] = (int) transformedExtent[0];
		dataExtent[3] = (int) transformedExtent[1];
		dataExtent[5] = (int) transformedExtent[2];

		for (idx = 0; idx < 6; idx += 2)
		{
			if (dataExtent[idx] > dataExtent[idx+1]) 
			{
				temp = dataExtent[idx];
				dataExtent[idx] = dataExtent[idx+1];
				dataExtent[idx+1] = temp;
			}
		}

		for (idx = 0; idx < 6; idx += 2)
		{
			// do the slide to 000 origin by subtracting the minimum extent
			inExtent[idx] += dataExtent[idx];
			inExtent[idx+1] += dataExtent[idx];
		}

		transformedExtent[0] = inExtent[0];
		transformedExtent[1] = inExtent[2];
		transformedExtent[2] = inExtent[4];
		this->Transform->GetLinearInverse()->TransformPoint(transformedExtent, 
			transformedExtent);
		outExtent[0] = (int) transformedExtent[0];
		outExtent[2] = (int) transformedExtent[1];
		outExtent[4] = (int) transformedExtent[2];

		transformedExtent[0] = inExtent[1];
		transformedExtent[1] = inExtent[3];
		transformedExtent[2] = inExtent[5];
		this->Transform->GetLinearInverse()->TransformPoint(transformedExtent,
			transformedExtent);
		outExtent[1] = (int) transformedExtent[0];
		outExtent[3] = (int) transformedExtent[1];
		outExtent[5] = (int) transformedExtent[2];

		for (idx = 0; idx < 6; idx += 2)
		{
			if (outExtent[idx] > outExtent[idx+1]) 
			{
				temp = outExtent[idx];
				outExtent[idx] = outExtent[idx+1];
				outExtent[idx+1] = temp;
			}
		}
	}

	vtkDebugMacro(<< "Inverse Transformed extent are:" 
		<< outExtent[0] << ", " << outExtent[1] << ", "
		<< outExtent[2] << ", " << outExtent[3] << ", "
		<< outExtent[4] << ", " << outExtent[5]);
}

void vtkOOCMetaImageReader::ComputeTransformedIncrements(vtkIdType inIncr[3],
														 vtkIdType outIncr[3])
{
	double transformedIncr[3];

	if (!this->Transform)
	{
		memcpy (outIncr, inIncr, 3 * sizeof (vtkIdType));
	}
	else
	{
		transformedIncr[0] = inIncr[0];
		transformedIncr[1] = inIncr[1];
		transformedIncr[2] = inIncr[2];
		this->Transform->TransformVector(transformedIncr, transformedIncr);
		outIncr[0] = (vtkIdType) transformedIncr[0];
		outIncr[1] = (vtkIdType) transformedIncr[1];
		outIncr[2] = (vtkIdType) transformedIncr[2];
		vtkDebugMacro(<< "Transformed Incr are:" 
			<< outIncr[0] << ", " << outIncr[1] << ", " << outIncr[2]);
	}
}


void vtkOOCMetaImageReader::ComputeInverseTransformedIncrements(vtkIdType inIncr[3],
																vtkIdType outIncr[3])
{
	double transformedIncr[3];

	if (!this->Transform)
	{
		memcpy (outIncr, inIncr, 3 * sizeof (vtkIdType));
	}
	else
	{
		transformedIncr[0] = inIncr[0];
		transformedIncr[1] = inIncr[1];
		transformedIncr[2] = inIncr[2];
		this->Transform->GetLinearInverse()->TransformVector(transformedIncr,
			transformedIncr);
		outIncr[0] = (vtkIdType) transformedIncr[0];
		outIncr[1] = (vtkIdType) transformedIncr[1];
		outIncr[2] = (vtkIdType) transformedIncr[2];
		vtkDebugMacro(<< "Inverse Transformed Incr are:" 
			<< outIncr[0] << ", " << outIncr[1] << ", " << outIncr[2]);
	}
}

//----------------------------------------------------------------------------
int vtkOOCMetaImageReader::CanReadFile(const char* fname)
{
        std::string filename = fname;
	if( filename == "" )
	{
		return false;
	}

	bool extensionFound = false;
        std::string::size_type mhaPos = filename.rfind(".mha");
        if ((mhaPos != std::string::npos)
		&& (mhaPos == filename.length() - 4))
	{
		extensionFound = true;
	}
        std::string::size_type mhdPos = filename.rfind(".mhd");
        if ((mhdPos != std::string::npos)
		&& (mhdPos == filename.length() - 4))
	{
		extensionFound = true;
	}
	if( !extensionFound )
	{
		return false;
	}

	// Now check the file content
	std::ifstream inputStream;

#ifdef __sgi
	inputStream.open( fname, ios::in );
#else
	inputStream.open( fname, ios::in | ios::binary );
#endif

	if( inputStream.fail() )
	{
		return false;
	}

	char key[8000];

	inputStream >> key;

	if( inputStream.eof() )
	{
		inputStream.close();
		return false;
	}

	if( strcmp(key,"NDims")==0 ) 
	{
		inputStream.close();
		return 3;
	}
	if( strcmp(key,"ObjectType")==0 ) 
	{
		inputStream.close();
		return 3;
	}
	if( strcmp(key,"TransformType")==0 ) 
	{
		inputStream.close();
		return 3;
	}
	if( strcmp(key,"ID")==0 ) 
	{
		inputStream.close();
		return 3;
	}
	if( strcmp(key,"ParentID")==0 ) 
	{
		inputStream.close();
		return 3;
	}
	if( strcmp(key,"BinaryData")==0 ) 
	{
		inputStream.close();
		return 3;
	}
	if( strcmp(key,"Comment")==0 ) 
	{
		inputStream.close();
		return 3;
	}
	if( strcmp(key,"AcquisitionDate")==0 ) 
	{
		inputStream.close();
		return 3;
	}
	if( strcmp(key,"Modality")==0 ) 
	{
		inputStream.close();
		return 3;
	}

	inputStream.close();
	return false;
}

//----------------------------------------------------------------------------
int vtkOOCMetaImageReader::GetDataByteOrder(void)
{
	return vtkmetaio::MET_SystemByteOrderMSB();
}
