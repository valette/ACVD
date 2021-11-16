include_directories(
    ${VTKSURFACE_INCLUDE_DIR}
)

option(USE_LLOYDCLUSTERING "Use Lloyd relaxations instead of new algorithm for clustering" OFF)
mark_as_advanced(USE_LLOYDCLUSTERING)


if(USE_LLOYDCLUSTERING)
    add_definitions (-D DOLloydClustering)
endif(USE_LLOYDCLUSTERING)

file( GLOB _vtkDiscreteRemeshing_hdrs *.h )
install( FILES ${_vtkDiscreteRemeshing_hdrs}
  DESTINATION ${INSTALL_INCLUDE_DIR}/ACVD/DiscreteRemeshing
)

if( BUILD_EXAMPLES )
  add_subdirectory( Examples )
endif()

## eof - CMakeLists.txt
