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
add_library(
vtkDiscreteRemeshing
vtkTestUniformClustering.cxx
xmlParser.cpp
vtkManifoldSimplification.cxx
)

target_link_libraries(vtkDiscreteRemeshing vtkSurface)

if( BUILD_EXAMPLES )
  add_subdirectory( Examples )
endif()

## eof - CMakeLists.txt
