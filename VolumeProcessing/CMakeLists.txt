add_library(vtkVolumeProcessing
vtkImageDataCleanLabels.cxx
)

target_link_libraries(vtkVolumeProcessing ${VTK_LIBRARIES})

if( BUILD_EXAMPLES )
  add_subdirectory( Examples )
endif()
