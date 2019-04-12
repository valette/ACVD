# Public headers
file(GLOB _vtkVolumeProcessing_hdrs include/*.h)

add_library(vtkVolumeProcessing
    src/vtkImageDataCleanLabels.cxx
    src/vtkOOCMetaImageReader.cxx
)

target_include_directories(vtkVolumeProcessing
    PUBLIC
      $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
      $<INSTALL_INTERFACE:${INSTALL_INCLUDE_DIR}/ACVD/VolumeProcessing>
)

set_target_properties(vtkVolumeProcessing
  PROPERTIES
    DEBUG_POSTFIX "-d"
    VERSION "${ACVD_VERSION}"
    SOVERSION "${ACVD_LIB_SOVERSION}"
    PUBLIC_HEADER "${_vtkVolumeProcessing_hdrs}"
)

install(
  TARGETS vtkVolumeProcessing
  EXPORT ACVDTargets
  RUNTIME DESTINATION ${INSTALL_BIN_DIR}
  LIBRARY DESTINATION ${INSTALL_LIB_DIR}
  ARCHIVE DESTINATION ${INSTALL_LIB_DIR}
  PUBLIC_HEADER DESTINATION ${INSTALL_INCLUDE_DIR}/ACVD/VolumeProcessing
)

if( BUILD_EXAMPLES )
  add_subdirectory( Examples )
endif()
