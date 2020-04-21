include_directories( ${ACVD_SOURCE_DIR}/Common )

set(COMMON_EXAMPLES
  viewer
  viewer2
  RandomTriangulation
  ExamplevtkSurface
  CleanMesh
  mesh2vtk
  vtk2ply
  stl2ply
  CheckOrientation
)

foreach(loop_var ${COMMON_EXAMPLES})

  if (VTK_VERSION VERSION_LESS "8.90.0")
    # old system
    include(${VTK_USE_FILE})
    add_executable(${loop_var} ${loop_var}.cxx)
    target_link_libraries(${loop_var} vtkSurface  ${VTK_LIBRARIES})
   else ()
    add_executable(${loop_var} ${loop_var}.cxx)
    target_link_libraries(${loop_var} vtkSurface ${VTK_LIBRARIES})
    vtk_module_autoinit(
    TARGETS ${loop_var}
    MODULES ${VTK_LIBRARIES}
  )
  endif () 

endforeach(loop_var)


install( TARGETS ${COMMON_EXAMPLES}
  EXPORT ACVDTargets
  RUNTIME DESTINATION ${INSTALL_BIN_DIR}
  LIBRARY DESTINATION ${INSTALL_LIB_DIR}
  ARCHIVE DESTINATION ${INSTALL_LIB_DIR}
)

