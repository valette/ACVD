set(lib_sources
	vtkImageDataCleanLabels.cxx
	vtkOOCMetaImageReader.cxx
)

if (VTK_VERSION VERSION_LESS "8.90.0")
  # old system
  include(${VTK_USE_FILE})
  add_library(vtkVolumeProcessing ${lib_sources} )
  target_link_libraries(vtkVolumeProcessing PRIVATE ${VTK_LIBRARIES})
else ()
  # include all components
  add_library(vtkVolumeProcessing ${lib_sources} )
  target_link_libraries(vtkVolumeProcessing PRIVATE ${VTK_LIBRARIES})

  # vtk_module_autoinit is needed
  vtk_module_autoinit(
    TARGETS vtkVolumeProcessing
    MODULES ${VTK_LIBRARIES}
    )
endif () 


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
