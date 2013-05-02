ADD_LIBRARY(vtkVolumeProcessing
vtkImageDataCleanLabels.cxx
)

TARGET_LINK_LIBRARIES(vtkVolumeProcessing vtkHybrid)

#-----------------------------------------------------------------------------
# Include the file 'InstallLibrary.cmake' which does the clean install job
SET(LIBRARY_NAME vtkVolumeProcessing)
SET(LIBRARY_NAME_BUILD_TREE_CONFIG ${VTKSURFACE_BINARY_DIR})
# Installation
# Headers
INSTALL_FILES(/include/${LIBRARY_NAME} "\\.h$")
# Binary
INSTALL_TARGETS(/lib ${LIBRARY_NAME})
SET(INSTALL_LIBRARY_FILES_DIR ${VTKSURFACE_SOURCE_DIR}/CMake)
INCLUDE(../CMake/InstallLibrary.cmake)
#-----------------------------------------------------------------------------


ADD_EXECUTABLE(VolumeCleanLabels VolumeCleanLabels.cxx)
TARGET_LINK_LIBRARIES(VolumeCleanLabels vtkHybrid vtkVolumeProcessing)

ADD_EXECUTABLE(VolumeSlice VolumeSlice.cxx)
TARGET_LINK_LIBRARIES(VolumeSlice vtkHybrid)

ADD_EXECUTABLE(Minc2Mhd Minc2Mhd.cxx)
TARGET_LINK_LIBRARIES(Minc2Mhd vtkHybrid)

ADD_EXECUTABLE(VolumeSubsample VolumeSubsample.cxx)
TARGET_LINK_LIBRARIES(VolumeSubsample vtkHybrid)

ADD_EXECUTABLE(VolumeCrop VolumeCrop.cxx)
TARGET_LINK_LIBRARIES(VolumeCrop vtkHybrid)

ADD_EXECUTABLE(VolumeMedian VolumeMedian.cxx)
TARGET_LINK_LIBRARIES(VolumeMedian vtkHybrid)

ADD_EXECUTABLE(VolumeAnisotropicDiffusion VolumeAnisotropicDiffusion.cxx)
TARGET_LINK_LIBRARIES(VolumeAnisotropicDiffusion vtkHybrid)

#ADD_EXECUTABLE(OOCRead OOCRead.cxx vtkOOCMetaImageReader.cxx BinarySubimageWriter.cxx ResampledVolumeGenerator.cxx vtkMHDFileWriter.cxx)
#TARGET_LINK_LIBRARIES(OOCRead vtkHybrid)

ADD_EXECUTABLE(png2raw png2raw.cxx)
TARGET_LINK_LIBRARIES(png2raw vtkHybrid)

ADD_EXECUTABLE(readimage ReadImage.cxx)
TARGET_LINK_LIBRARIES(readimage vtkHybrid)
