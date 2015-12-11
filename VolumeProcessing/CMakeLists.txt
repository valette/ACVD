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

if( BUILD_EXAMPLES )
	add_subdirectory( Examples )
endif()
