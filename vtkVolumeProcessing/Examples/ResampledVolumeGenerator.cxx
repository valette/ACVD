#include "ResampledVolumeGenerator.h"

ResampledVolumeGenerator::ResampledVolumeGenerator ( ) {
	metaReader = vtkOOCMetaImageReader::New ( );
	resampler = vtkImageResample::New ( );

	_divisionFactorX = 0;
	_divisionFactorY = 0;
	_divisionFactorZ = 0;

	_resamplingFactorX = 0;
	_resamplingFactorY = 0;
	_resamplingFactorZ = 0;

	_dimensionX = 0;
	_dimensionY = 0;
	_dimensionZ = 0;

	_stepX = 0;
	_stepY = 0;
	_stepZ = 0;

	_dataType = 0;
	_dataTypeSize = 0;

	_write = true;

	//writer = new BinarySubimageWriter();
	//mhdwriter = new vtkMHDFileWriter();

}

ResampledVolumeGenerator::~ResampledVolumeGenerator ( ) {
	metaReader->Delete();
	resampler->Delete();
}

void ResampledVolumeGenerator::SetDivisionFactors ( int divisionFactorX, int divisionFactorY, int divisionFactorZ ) {
	_divisionFactorX = divisionFactorX;
	_divisionFactorY = divisionFactorY;
	_divisionFactorZ = divisionFactorZ;
}

void ResampledVolumeGenerator::SetResamplingFactors ( double resamplingFactorX, double resamplingFactorY, double resamplingFactorZ ){
	_resamplingFactorX = resamplingFactorX;
	_resamplingFactorY = resamplingFactorY;
	_resamplingFactorZ = resamplingFactorZ;
}

void ResampledVolumeGenerator::SetFileName( char *filename ){

	strcpy(_filename, filename);
	strcpy(_originalFilename, _filename);
	metaReader->SetFileName(_originalFilename);

	//Truncate the filename as to eliminate .mhd
	ObtainTruncatedFilename( );

	//Add the correct ending to the filename
	RawDataFilename();

	//std::cout << "Filename that's to be assigned " << _filename << std::endl;
	//Assign the filename to the writer

	//Get the element data file out of the modified filename
	AssignElementDataFile();
}

void ResampledVolumeGenerator::AssignElementDataFile ( ){

	strcpy(_elementDataName, _rawFilename );
	//char* pointer = strrchr(_filename, '/');
	//strcpy(_elementDataName, pointer);
	//std::cout << "Delete this shit, element data name " << _elementDataName << std::endl;
}

void ResampledVolumeGenerator::ObtainTruncatedFilename ( ){
	int i;
	long l = static_cast<long>( strlen( _filename ) );
	_truncatedFilename = new char[l+4];

	for(i=l; i>=0; i--){
		if( _filename[i] == '.' ){
			strcpy(_truncatedFilename, _filename);
			_truncatedFilename[i] = '\0';
			break;
		}
	}
}

void ResampledVolumeGenerator::RawDataFilename(){
	char *tempFilename;
	int strLn = strlen( _truncatedFilename );
	tempFilename = new char[ strLn + 8 ];
	strcpy(tempFilename, _truncatedFilename );
	tempFilename[strLn] = '_';
	tempFilename[strLn+1] = 'a';
	tempFilename[strLn + 2] = 'l';
	tempFilename[strLn + 3] = 't';
	tempFilename[strLn + 4] = '.';
	tempFilename[strLn + 5] = 'r';
	tempFilename[strLn + 6] = 'a';
	tempFilename[strLn + 7] = 'w';
	tempFilename[strLn + 8] = '\0';
	strcpy(_rawFilename, tempFilename);
}

void ResampledVolumeGenerator::MHDFilename ( ) {
	char *tempFilename;
	int strLn = strlen( _truncatedFilename );
	tempFilename = new char[ strLn + 8 ];
	strcpy(tempFilename, _truncatedFilename );
	tempFilename[strLn] = '_';
	tempFilename[strLn+1] = 'a';
	tempFilename[strLn + 2] = 'l';
	tempFilename[strLn + 3] = 't';
	tempFilename[strLn + 4] = '.';
	tempFilename[strLn + 5] = 'm';
	tempFilename[strLn + 6] = 'h';
	tempFilename[strLn + 7] = 'd';
	tempFilename[strLn + 8] = '\0';
	strcpy(_mhdFilename, tempFilename);
}

void ResampledVolumeGenerator::GetImageInformation ( ){
	metaReader->ReadImageParameters ( );
	int dimensions[6];
	metaReader->GetDataExtent(dimensions);

	//Dimensions = extent + 1
	_dimensionX = dimensions[1]+1;
	_dimensionY = dimensions[3]+1;
	_dimensionZ = dimensions[5]+1;

	//Calculate the dimensions of the resampled image
	_resampledDimensionX = floor ( (double) ( (_dimensionX) * _resamplingFactorX) );
	_resampledDimensionY = floor ( (double) ( (_dimensionY) * _resamplingFactorY) );
	_resampledDimensionZ = floor ( (double) ( (_dimensionZ) * _resamplingFactorZ) );

	//Calculate the subvolume step that's to be used
	_stepX = ceil((double)_dimensionX/_divisionFactorX);
	_stepY = ceil((double)_dimensionY/_divisionFactorY);
	_stepZ = ceil((double)_dimensionZ/_divisionFactorZ);


	//Obtain the data type and its respective size
	_dataType = metaReader->GetDataScalarType ( );
	int numberOfChannels = metaReader->GetNumberOfScalarComponents();

	//Create a vector with the resampled dimensions
	int _rDimensions[3];
	_rDimensions[0] = _resampledDimensionX; 
	_rDimensions[1] = _resampledDimensionY;
	_rDimensions[2] = _resampledDimensionZ;

	//Set the information required for correctly creating the MHD File
	double spacing[3];
	double position[3];
	//int *centerOfRotation;
	//int *offset;

	metaReader->GetDataSpacing(spacing);
	metaReader->GetDataOrigin(position);

	inter = vtkImageData::New();
	inter->SetOrigin(position);
	inter->SetSpacing(spacing);


	if( _write ){
		//Set the information required for correctly creating the raw data file
		writer->SetImageDimensions(_resampledDimensionX, _resampledDimensionY, _resampledDimensionZ );
		writer->SetPixelSizeAndType( _dataType );
		writer->SetNumberOfChannels( numberOfChannels );


		//Assign all the relevant data to the mhd file writer
		MHDFilename();
		mhdwriter->SetMHDFileName(_mhdFilename);
		mhdwriter->SetObjectType("Image");
		mhdwriter->SetNDimensions(metaReader->GetFileDimensionality());
		mhdwriter->SetDimensions(_resampledDimensionX, _resampledDimensionY, _resampledDimensionZ);
		mhdwriter->SetElementType(_dataType);
		mhdwriter->SetHeaderSize(0);
		mhdwriter->SetNumberOfChannels(numberOfChannels);
		mhdwriter->SetElementSpacing(spacing[0], spacing[1], spacing[2]);
		mhdwriter->SetPosition(position[0], position[1], position[2]);

		//Verifier avec Sébastien
		mhdwriter->SetCenterOfRotation(0,0,0);
		mhdwriter->SetOffset(0,0,0);

		mhdwriter->SetIsBinaryData(true);

		int orderBytes = metaReader->GetDataByteOrder();
		bool isMSB = true;
		if (orderBytes!=1)
			isMSB = false;
		mhdwriter->SetIsByteOrderMSB(isMSB);
		mhdwriter->SetIsCompressedData(false);
		mhdwriter->SetElementDataFile(_elementDataName);
	}
}
/*
void ResampledVolumeGenerator::ConfigureFilenames ( ) {
	strcpy(_filename, GetData());
	std::cout << "Filename " <<  _filename << std::endl;
	strcpy(_originalFilename, _filename);
	metaReader->SetFileName(_originalFilename);

	//Truncate the filename as to eliminate .mhd
	ObtainTruncatedFilename( );

	//Add the correct ending to the filename
	RawDataFilename();

	//std::cout << "Filename that's to be assigned " << _filename << std::endl;
	//Assign the filename to the writer

	//Get the element data file out of the modified filename
	AssignElementDataFile();
}
*/
void ResampledVolumeGenerator::StartWork ( ){
//	ConfigureFilenames ( );
	//if ( write ){
	std::cout << "Write" << std::endl;
	writer = new BinarySubimageWriter();
	mhdwriter = new vtkMHDFileWriter();
	writer->AssignFilename( _rawFilename );
	writer->OpenFile();
	//}
	GetImageInformation();
	//X
	int _actuX = 0; 
	int _finX = _actuX + _stepX;
	if( _finX >= _dimensionX )
		_finX = _dimensionX - 1;

	//Y
	int _actuY = 0; 
	int _finY = _actuY + _stepY;
	if ( _finY >= _dimensionY )
		_finY = _dimensionY - 1;

	//Z
	int _actuZ = 0;
	int _finZ = _actuZ + _stepZ;
	if( _finZ >= _dimensionZ )
		_finZ = _dimensionZ - 1;


	std::cout << "We're in for " << _divisionFactorX*_divisionFactorY*_divisionFactorZ << " iterations" << std::endl;
	std::cout << "Begin subdividing and resampling" << std::endl;
	int iteration = 1;
	for (int z = 0; z < _divisionFactorZ; z++) {
		for (int x = 0; x < _divisionFactorX; x++) {
			for( int y = 0; y < _divisionFactorY; y++ ){

				std::cout << "	Iteration number " << iteration << std::endl;
				iteration++;
				//ENCORE A FAIRE: Trouver une manière pour ne pas devoir créer un nouveau reader/image à chaque fois
				//ENCORE A FAIRE: Il y a quelque chose d'inefficient . J'utilise une instance du reader pour obtenir les paramètres de l'image,
				//et de nombreuses autres instances pour lire les subvolumes. Ce n'est pas très intelligent, est-ce qu'il n'existe pas une manière
				//plus intelligente de faire cela ?
				/*std::cout << "This is the filename we're working with " << _filename << std::endl;*/
				vtkOOCMetaImageReader *reader = vtkOOCMetaImageReader::New();
				vtkImageData* temp = vtkImageData::New();
				reader->SetFileName ( _originalFilename );

				int indexes [6] = {_actuX, _finX, _actuY, _finY, _actuZ, _finZ};
				std::cout << "		Indexes " << _actuX << " " << _finX << " " << _actuY << " " << _finY << " " << _actuZ << " " << _finZ << std::endl;
				reader->SetDataVOI (indexes);
				reader->Update();
				temp->ShallowCopy(reader->GetOutput());

				//ENCORE A FAIRE: Trouver une manière pour ne pas devoir créer un nouveau resampler/image à chaque fois
				//Voilà, quelqu'un d'autre est en charge du sous-échantillonnage!
				vtkImageResample *resampler = vtkImageResample::New();
				resampler->SetInput(temp);
				resampler->SetAxisMagnificationFactor(0, _resamplingFactorX);
				resampler->SetAxisMagnificationFactor(1, _resamplingFactorY);
				resampler->SetAxisMagnificationFactor(2, _resamplingFactorZ);
				resampler->Update();
				vtkImageData* sousVolume = vtkImageData::New();
				sousVolume->ShallowCopy(resampler->GetOutput());
				//imageMatrix[x][y][z] = *sousVolume;



				//ENCORE A FAIRE: Comment faire pour utiliser la même instance du reader pour obtenir la metadata et les données du fichier?
				//ENCORE A FAIRE: écrire ligne par ligne et pas pixel par pixel
				//En fait, tout ce qui vient après devrait se faire dans une autre classe que soit en charge 
				//de l'écriture dans le fichier!
				std::cout << "Is this trick writing?" << std::endl;
				int _sextent[6];
				sousVolume->GetExtent(_sextent);
				writer->ChangeSubimage(sousVolume);
				writer->SetSubimageExtent(_sextent);
				writer->WriteToFile();

				//ENCORE A FAIRE: Il faut bien regarder où l'on est en train de demander de la mémoire, ici et ailleurs. Effectivement,
				//la création constante des nouvelles classes fait péter le truc. C'est normal, les images de Sébastien sont géantes.
				//Du coup, il faut faire gaffe.
				reader->Delete();
				resampler->Delete();
				sousVolume->Delete();
				temp->Delete();

				if( _finY==(_dimensionY-1) ){
					_actuY = 0;
					_finY = _stepY;
				}else if ( (_finY+_stepY)>(_dimensionY-1) ){
					_actuY = _actuY + _stepY;
					_finY = (_dimensionY-1);
				}else{
					_actuY = _actuY + _stepY;
					_finY = _finY + _stepY;
				}
			}
			if( _finX==(_dimensionX-1) ){
				_actuX = 0;
				_finX = _stepX;
			}else if ( (_finX+_stepX)>(_dimensionX-1) ){
				_actuX = _actuX + _stepX;
				_finX = (_dimensionX-1);
			}else{
				_actuX = _actuX + _stepX;
				_finX = _finX + _stepX;
			}
		}
		if( _finZ==(_dimensionZ-1) ){
			break;
		}else if ( (_finZ+_stepZ)>(_dimensionZ-1) ){
			_actuZ = _actuZ + _stepZ;
			_finZ = (_dimensionZ-1);
		}else{
			_actuZ = _actuZ + _stepZ;
			_finZ = _finZ + _stepZ;
		}
	}

	//if( write ){
		std::cout << "End subdividing and resampling" << std::endl;
		writer->CloseFile();
		std::cout << "Begin writing the mhdFile" << std::endl;
		mhdwriter->WriteToFile();
	//}
	/*NotifyFinished();
	StopExecution();*/
}

void ResampledVolumeGenerator::UnthreadedExecute ( ){
	//if ( write ){
	//	writer = new BinarySubimageWriter();
	//	mhdwriter = new vtkMHDFileWriter();
	//	writer->AssignFilename( _rawFilename );
	//}
	StartWork();
}

char* ResampledVolumeGenerator::GetNewMHDFileName ( ) {
	return _mhdFilename;
}

void  ResampledVolumeGenerator::SetWriteToRawFile ( bool write ){
	_write = write ;
	std::cout << "Write " << _write << std::endl;
}
