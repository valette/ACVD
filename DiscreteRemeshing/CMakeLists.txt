include_directories(
    ${VTKSURFACE_INCLUDE_DIR}
)

option(USE_MULTITHREADING "Use Multithreaded classes" OFF)
mark_as_advanced(USE_MULTITHREADING)

option(USE_LLOYDCLUSTERING "Use Lloyd relaxations instead of new algorithm for clustering" OFF)
mark_as_advanced(USE_LLOYDCLUSTERING)


# --------------------------------------------------------------------
# Switch to enable Multithreaded clustering

if(USE_MULTITHREADING)
    add_definitions (-D DOmultithread)
endif(USE_MULTITHREADING)
if(USE_LLOYDCLUSTERING)
    add_definitions (-D DOLloydClustering)
endif(USE_LLOYDCLUSTERING)


 # --------------------------------------------------------------------------
 # Library compilation
set(lib_sources
  vtkTestUniformClustering.cxx
  vtkManifoldSimplification.cxx
)


if (VTK_VERSION VERSION_LESS "8.90.0")
  # old system
  include(${VTK_USE_FILE})
  add_library(vtkDiscreteRemeshing ${lib_sources} )
  target_link_libraries(vtkDiscreteRemeshing vtkSurface  ${VTK_LIBRARIES})
else ()
  # include all components
  add_library(vtkDiscreteRemeshing ${lib_sources} )
  target_link_libraries(vtkDiscreteRemeshing vtkSurface ${VTK_LIBRARIES})

  # vtk_module_autoinit is needed
  vtk_module_autoinit(
    TARGETS vtkDiscreteRemeshing
    MODULES ${VTK_LIBRARIES}
    )
endif () 


set_target_properties(vtkDiscreteRemeshing PROPERTIES
  DEBUG_POSTFIX "-d"
  VERSION "${ACVD_LIB_VERSION}"
  SOVERSION "${ACVD_LIB_SOVERSION}"
)

install( TARGETS vtkDiscreteRemeshing
  EXPORT ACVDTargets
  RUNTIME DESTINATION ${INSTALL_BIN_DIR}
  LIBRARY DESTINATION ${INSTALL_LIB_DIR}
  ARCHIVE DESTINATION ${INSTALL_LIB_DIR}
)

file( GLOB _vtkDiscreteRemeshing_hdrs *.h )
install( FILES ${_vtkDiscreteRemeshing_hdrs}
  DESTINATION ${INSTALL_INCLUDE_DIR}/ACVD/DiscreteRemeshing
)

if( BUILD_EXAMPLES )
  add_subdirectory( Examples )
endif()

## eof - CMakeLists.txt
