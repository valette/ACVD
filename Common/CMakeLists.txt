 
 # --------------------------------------------------------------------------
 # Library compilation
 ADD_LIBRARY(vtkSurface
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
 )
 

# Set the added libraries
SET(LIB_ADDED 
vtkCommon
vtkFiltering
vtkImaging
vtkIO
vtkRendering
vtkHybrid)

TARGET_LINK_LIBRARIES(vtkSurface ${LIB_ADDED})

SET(COMMON_EXAMPLES
viewer
viewer2
RandomTriangulation
ExamplevtkSurface
)
FOREACH(loop_var ${COMMON_EXAMPLES})
	ADD_EXECUTABLE(${loop_var} ${loop_var}.cxx)
	TARGET_LINK_LIBRARIES(${loop_var} vtkSurface
	vtkCommon vtkFiltering vtkGraphics vtkImaging vtkIO vtkRendering vtkHybrid)
ENDFOREACH(loop_var)

#-----------------------------------------------------------------------------
# Include the file 'InstallLibrary.cmake' which does the clean install job
SET(LIBRARY_NAME vtkSurface)
SET(LIBRARY_NAME_BUILD_TREE_CONFIG ${CMAKE_BINARY_DIR})
# Installation
# Headers
INSTALL_FILES(/include/${LIBRARY_NAME} "\\.h$")
# Binary
INSTALL_TARGETS(/lib ${LIBRARY_NAME})
SET(INSTALL_LIBRARY_FILES_DIR ${CMAKE_SOURCE_DIR}/CMake)
INCLUDE(../CMake/InstallLibrary.cmake)
#-----------------------------------------------------------------------------
