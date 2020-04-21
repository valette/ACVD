include_directories(
    ${VTK_INCLUDE_DIR}
    ${VTKSURFACE_INCLUDE_DIR}
    ${VTKDISCRETEREMESHING_INCLUDE_DIR}
)

set(DISCRETEREMESHING_EXAMPLES
ACVD
ACVDQ
AnisotropicRemeshing
AnisotropicRemeshingQ
VolumeAnalysis
ManifoldSimplification
)

foreach(loop_var ${DISCRETEREMESHING_EXAMPLES})
  if (VTK_VERSION VERSION_LESS "8.90.0")
    # old system
    include(${VTK_USE_FILE})
    add_executable(${loop_var} ${loop_var}.cxx)
    target_link_libraries(${loop_var} vtkDiscreteRemeshing  ${VTK_LIBRARIES})
  else ()
    add_executable(${loop_var} ${loop_var}.cxx)
    target_link_libraries(${loop_var} vtkDiscreteRemeshing ${VTK_LIBRARIES})
    vtk_module_autoinit(
    TARGETS ${loop_var}
    MODULES ${VTK_LIBRARIES}
  )
  endif () 
endforeach(loop_var)

install( TARGETS ${DISCRETEREMESHING_EXAMPLES}
  EXPORT ACVDTargets
  RUNTIME DESTINATION ${INSTALL_BIN_DIR}
  LIBRARY DESTINATION ${INSTALL_LIB_DIR}
  ARCHIVE DESTINATION ${INSTALL_LIB_DIR}
)
