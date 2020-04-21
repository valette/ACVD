include_directories( ${ACVD_SOURCE_DIR}/VolumeProcessing )

set(VolumeProcessing_EXAMPLES
	img2raw
	Minc2Mhd
	png2raw
	ReadImage
	VolumeAnisotropicDiffusion
	VolumeCleanLabels
	VolumeCrop
	VolumeMedian
	VolumeOOCSlice
	VolumeSlice
	VolumeSubsample
)

foreach(loop_var ${VolumeProcessing_EXAMPLES})
  if (VTK_VERSION VERSION_LESS "8.90.0")
    # old system
    include(${VTK_USE_FILE})
    add_executable(${loop_var} ${loop_var}.cxx)
    target_link_libraries(${loop_var} vtkSurface vtkVolumeProcessing  ${VTK_LIBRARIES})
  else ()
    add_executable(${loop_var} ${loop_var}.cxx)
    target_link_libraries(${loop_var} vtkSurface vtkVolumeProcessing ${VTK_LIBRARIES})
    vtk_module_autoinit(
    TARGETS ${loop_var}
    MODULES ${VTK_LIBRARIES}
  )
  endif () 
endforeach(loop_var)

install( TARGETS ${VolumeProcessing_EXAMPLES}
  EXPORT ACVDTargets
  RUNTIME DESTINATION ${INSTALL_BIN_DIR}
  LIBRARY DESTINATION ${INSTALL_LIB_DIR}
  ARCHIVE DESTINATION ${INSTALL_LIB_DIR}
)
