ADD_LIBRARY(vtkVolumeProcessing
vtkImageDataCleanLabels.cxx
)

TARGET_LINK_LIBRARIES(vtkVolumeProcessing ${VTK_LIBRARIES})

ADD_EXECUTABLE(VolumeCleanLabels VolumeCleanLabels.cxx)
TARGET_LINK_LIBRARIES(VolumeCleanLabels vtkVolumeProcessing ${VTK_LIBRARIES})

ADD_EXECUTABLE(VolumeSlice VolumeSlice.cxx)
TARGET_LINK_LIBRARIES(VolumeSlice ${VTK_LIBRARIES})

ADD_EXECUTABLE(Minc2Mhd Minc2Mhd.cxx)
TARGET_LINK_LIBRARIES(Minc2Mhd ${VTK_LIBRARIES})

ADD_EXECUTABLE(VolumeSubsample VolumeSubsample.cxx)
TARGET_LINK_LIBRARIES(VolumeSubsample ${VTK_LIBRARIES})

ADD_EXECUTABLE(VolumeCrop VolumeCrop.cxx)
TARGET_LINK_LIBRARIES(VolumeCrop ${VTK_LIBRARIES})

ADD_EXECUTABLE(VolumeMedian VolumeMedian.cxx)
TARGET_LINK_LIBRARIES(VolumeMedian ${VTK_LIBRARIES})

ADD_EXECUTABLE(VolumeAnisotropicDiffusion VolumeAnisotropicDiffusion.cxx)
TARGET_LINK_LIBRARIES(VolumeAnisotropicDiffusion ${VTK_LIBRARIES})

#ADD_EXECUTABLE(OOCRead OOCRead.cxx vtkOOCMetaImageReader.cxx BinarySubimageWriter.cxx ResampledVolumeGenerator.cxx vtkMHDFileWriter.cxx)
#TARGET_LINK_LIBRARIES(OOCRead vtkHybrid)

ADD_EXECUTABLE(png2raw png2raw.cxx)
TARGET_LINK_LIBRARIES(png2raw ${VTK_LIBRARIES})

ADD_EXECUTABLE(readimage ReadImage.cxx)
TARGET_LINK_LIBRARIES(readimage ${VTK_LIBRARIES})
