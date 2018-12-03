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
  add_executable(${loop_var} ${loop_var}.cxx)
  target_link_libraries(${loop_var} vtkDiscreteRemeshing )
endforeach(loop_var)

install( TARGETS ${DISCRETEREMESHING_EXAMPLES}
  RUNTIME DESTINATION ${INSTALL_BIN_DIR}
  LIBRARY DESTINATION ${INSTALL_LIB_DIR}
  ARCHIVE DESTINATION ${INSTALL_LIB_DIR}
)
