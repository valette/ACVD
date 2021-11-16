 # --------------------------------------------------------------------------
 # Library compilation
add_library(vtkSurface
 RenderWindow.cxx
 vtkDelaunay.cxx
 vtkDualMeshDisplay.cxx
 vtkManifoldSimplification.cxx
 vtkOFFReader.cxx
 vtkOFFWriter.cxx
 vtkQuadricTools.cxx
 vtkRandomTriangulation.cxx
 vtkSurface.cxx
 vtkSurfaceBase.cxx
 vtkVolumeProperties.cxx
 vtkNeighbourhoodComputation.cxx
 vtkCurvatureMeasure.cxx
 )
target_link_libraries(vtkSurface ${VTK_LIBRARIES})

set_target_properties( vtkSurface PROPERTIES
  DEBUG_POSTFIX "-d"
  VERSION "${ACVD_LIB_VERSION}"
  SOVERSION "${ACVD_LIB_SOVERSION}"
)

install( TARGETS vtkSurface
  EXPORT ACVDTargets
  RUNTIME DESTINATION ${INSTALL_BIN_DIR}
  LIBRARY DESTINATION ${INSTALL_LIB_DIR}
  ARCHIVE DESTINATION ${INSTALL_LIB_DIR}
)

file( GLOB _vtkSurface_hdrs *.h )

install( FILES ${_vtkSurface_hdrs}
  DESTINATION ${INSTALL_INCLUDE_DIR}/ACVD/Common
)

if( BUILD_EXAMPLES )
  add_subdirectory( Examples )
endif()
