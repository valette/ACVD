include_directories(
    ${VTK_INCLUDE_DIR}
    ${VTKSURFACE_INCLUDE_DIR}
    ${VTKDISCRETEREMESHING_INCLUDE_DIR}
)

set(DISCRETEREMESHING_EXAMPLES
ACVD
ACVDP
ACVDQ
ACVDQP
AnisotropicRemeshing
AnisotropicRemeshingQ
AnisotropicRemeshingQP
VolumeAnalysis
)

find_package(OpenMP)

if (OPENMP_FOUND)
set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif (OPENMP_FOUND)

foreach(loop_var ${DISCRETEREMESHING_EXAMPLES})
  add_executable(${loop_var} ${loop_var}.cxx)
  target_link_libraries(${loop_var} vtkSurface ${VTK_LIBRARIES} )
endforeach(loop_var)

install( TARGETS ${DISCRETEREMESHING_EXAMPLES}
  EXPORT ACVDTargets
  RUNTIME DESTINATION ${INSTALL_BIN_DIR}
  LIBRARY DESTINATION ${INSTALL_LIB_DIR}
  ARCHIVE DESTINATION ${INSTALL_LIB_DIR}
)
