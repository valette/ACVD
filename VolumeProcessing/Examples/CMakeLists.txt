include_directories( ${ACVD_SOURCE_DIR}/VolumeProcessing )

add_executable(VolumeCleanLabels VolumeCleanLabels.cxx)
target_link_libraries(VolumeCleanLabels vtkVolumeProcessing ${VTK_LIBRARIES})

add_executable(VolumeSlice VolumeSlice.cxx)
target_link_libraries(VolumeSlice ${VTK_LIBRARIES})

add_executable(VolumeOOCSlice VolumeOOCSlice.cxx)
target_link_libraries(VolumeOOCSlice vtkVolumeProcessing ${VTK_LIBRARIES})

add_executable(Minc2Mhd Minc2Mhd.cxx)
target_link_libraries(Minc2Mhd ${VTK_LIBRARIES})

add_executable(VolumeSubsample VolumeSubsample.cxx)
target_link_libraries(VolumeSubsample ${VTK_LIBRARIES})

add_executable(VolumeCrop VolumeCrop.cxx)
target_link_libraries(VolumeCrop ${VTK_LIBRARIES})

add_executable(VolumeMedian VolumeMedian.cxx)
target_link_libraries(VolumeMedian ${VTK_LIBRARIES})

add_executable(VolumeAnisotropicDiffusion VolumeAnisotropicDiffusion.cxx)
target_link_libraries(VolumeAnisotropicDiffusion ${VTK_LIBRARIES})

#add_executable(OOCRead OOCRead.cxx vtkOOCMetaImageReader.cxx BinarySubimageWriter.cxx ResampledVolumeGenerator.cxx vtkMHDFileWriter.cxx)
#target_link_libraries(OOCRead vtkHybrid)

add_executable(png2raw png2raw.cxx)
target_link_libraries(png2raw ${VTK_LIBRARIES})

add_executable(readimage ReadImage.cxx)
target_link_libraries(readimage ${VTK_LIBRARIES})

install( TARGETS VolumeCleanLabels VolumeSlice Minc2Mhd VolumeSubsample VolumeCrop VolumeMedian VolumeAnisotropicDiffusion png2raw readimage
  EXPORT ACVDTargets
  RUNTIME DESTINATION ${INSTALL_BIN_DIR}
  LIBRARY DESTINATION ${INSTALL_LIB_DIR}
  ARCHIVE DESTINATION ${INSTALL_LIB_DIR}
)
