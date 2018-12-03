# --------------------------------------------------------------------
# Switch to enable Multithreaded clustering

option(USE_MULTITHREADING "Use Multithreaded classes" OFF)
mark_as_advanced(USE_MULTITHREADING)

option(USE_LLOYDCLUSTERING "Use Lloyd relaxations instead of new algorithm for clustering" OFF)
mark_as_advanced(USE_LLOYDCLUSTERING)

if(USE_MULTITHREADING)
    add_compile_definitions(DOmultithread)
endif(USE_MULTITHREADING)
if(USE_LLOYDCLUSTERING)
    add_compile_definitions(DOLloydClustering)
endif(USE_LLOYDCLUSTERING)


# --------------------------------------------------------------------------
# Library compilation

# Public headers
file( GLOB _vtkDiscreteRemeshing_hdrs include/*.h )

# Setup library
add_library(vtkDiscreteRemeshing
    src/vtkTestUniformClustering.cxx
    src/vtkManifoldSimplification.cxx
)
target_include_directories(vtkDiscreteRemeshing
    PUBLIC
      $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
      $<INSTALL_INTERFACE:${INSTALL_INCLUDE_DIR}/ACVD/DiscreteRemeshing>
)
target_link_libraries(vtkDiscreteRemeshing
    PUBLIC
        vtkSurface
)
set_target_properties(vtkDiscreteRemeshing
  PROPERTIES
    DEBUG_POSTFIX "-d"
    VERSION "${ACVD_VERSION}"
    SOVERSION "${ACVD_LIB_SOVERSION}"
    PUBLIC_HEADER "${_vtkDiscreteRemeshing_hdrs}"
)

install(
  TARGETS vtkDiscreteRemeshing
  EXPORT ACVDTargets
  RUNTIME DESTINATION ${INSTALL_BIN_DIR}
  LIBRARY DESTINATION ${INSTALL_LIB_DIR}
  ARCHIVE DESTINATION ${INSTALL_LIB_DIR}
  PUBLIC_HEADER DESTINATION ${INSTALL_INCLUDE_DIR}/ACVD/DiscreteRemeshing
)

if( BUILD_EXAMPLES )
  add_subdirectory( Examples )
endif()
