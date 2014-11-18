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
  add_executable(${loop_var} ${loop_var}.cxx)
  target_link_libraries(${loop_var} vtkSurface ${VTK_LIBRARIES})
endforeach(loop_var)

