add_library(vtkVolumeProcessing
vtkImageDataCleanLabels.cxx
vtkOOCMetaImageReader.cxx
)

target_link_libraries(vtkVolumeProcessing ${VTK_LIBRARIES})


set_target_properties( vtkVolumeProcessing PROPERTIES
  DEBUG_POSTFIX "-d"
  VERSION "${ACVD_LIB_VERSION}"
  SOVERSION "${ACVD_LIB_SOVERSION}"
)

install( TARGETS vtkVolumeProcessing
  EXPORT ACVDTargets
  RUNTIME DESTINATION ${INSTALL_BIN_DIR}
  LIBRARY DESTINATION ${INSTALL_LIB_DIR}
  ARCHIVE DESTINATION ${INSTALL_LIB_DIR}
)

file( GLOB _vtkVolumeProcessing_hdrs *.h )

install( FILES ${_vtkVolumeProcessing_hdrs}
  DESTINATION ${INSTALL_INCLUDE_DIR}/ACVD/VolumeProcessing
)

<<<<<<< HEAD
ADD_EXECUTABLE(readimage ReadImage.cxx)
TARGET_LINK_LIBRARIES(readimage vtkHybrid)

ADD_EXECUTABLE(VolumeOOCSlice VolumeOOCSlice.cxx)
TARGET_LINK_LIBRARIES(VolumeOOCSlice vtkHybrid vtkVolumeProcessing)

ADD_EXECUTABLE(img2raw img2raw.cxx)
TARGET_LINK_LIBRARIES(img2raw vtkHybrid)
=======
if( BUILD_EXAMPLES )
  add_subdirectory( Examples )
endif()
>>>>>>> d6582bd360bea6b52c1fe32a32c2509d03fba282
