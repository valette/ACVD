 # --------------------------------------------------------------------------
 # Library compilation
add_library(vtkSurface
 RenderWindow.cxx
 vtkDelaunay.cxx
 vtkDualMeshDisplay.cxx
 vtkOFFReader.cxx
 vtkOFFWriter.cxx
 vtkQuadricTools.cxx
 vtkRandomTriangulation.cxx
 vtkSMFReader.cxx
 vtkSurface.cxx
 vtkSurfaceBase.cxx
 vtkVolumeProperties.cxx
 vtkNeighbourhoodComputation.cxx
 vtkCurvatureMeasure.cxx
 vtkMyMinimalStandardRandomSequence.cxx
#  vtkJSWriter.cxx
 )
target_link_libraries(vtkSurface ${VTK_LIBRARIES})

if( BUILD_EXAMPLES )
  add_subdirectory( Examples )
endif()
