 # --------------------------------------------------------------------------
 # Library compilation

set(lib_sources
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

if (VTK_VERSION VERSION_LESS "8.90.0")
  # old system
  include(${VTK_USE_FILE})
  add_library(vtkSurface ${lib_sources} )
  target_link_libraries(vtkSurface ${VTK_LIBRARIES})
else ()
  # include all components
  add_library(vtkSurface ${lib_sources} )
  target_link_libraries(vtkSurface ${VTK_LIBRARIES})

  # vtk_module_autoinit is needed
  vtk_module_autoinit(
    TARGETS vtkSurface
    MODULES ${VTK_LIBRARIES}
    )
endif () 


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
