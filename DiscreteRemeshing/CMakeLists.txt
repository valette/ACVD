INCLUDE_DIRECTORIES(
    ${VTKSURFACE_INCLUDE_DIR}
)
 
 
 # --------------------------------------------------------------------------
 # Library compilation
 ADD_LIBRARY(
vtkDiscreteRemeshing 
vtkTestUniformClustering.cxx
)

SET(LIB_ADDED
vtkSurface
vtkCommon
vtkFiltering
vtkImaging
vtkIO
vtkRendering
vtkHybrid)

TARGET_LINK_LIBRARIES(vtkDiscreteRemeshing ${LIB_ADDED})

#-----------------------------------------------------------------------------
# Include the file 'InstallLibrary.cmake' which does the clean install job
SET(LIBRARY_NAME vtkDiscreteRemeshing)
SET(LIBRARY_NAME_BUILD_TREE_CONFIG ${CMAKE_BINARY_DIR})
# Installation
# Headers
INSTALL_FILES(/include/${LIBRARY_NAME} "\\.h$")
# Binary
INSTALL_TARGETS(/lib ${LIBRARY_NAME})
SET(INSTALL_LIBRARY_FILES_DIR ${CMAKE_SOURCE_DIR}/CMake)
INCLUDE(../CMake/InstallLibrary.cmake)
#-----------------------------------------------------------------------------

INCLUDE_DIRECTORIES(
    ${VTK_INCLUDE_DIR}
    ${VTKSURFACE_INCLUDE_DIR}
    ${VTKDISCRETEREMESHING_INCLUDE_DIR}
)

SET(DISCRETEREMESHING_EXAMPLES
ACVD
ACVDQ
AnisotropicRemeshingQ
ManifoldSimplification
)

FOREACH(loop_var ${DISCRETEREMESHING_EXAMPLES})
	ADD_EXECUTABLE(${loop_var} ${loop_var}.cxx)
	TARGET_LINK_LIBRARIES(${loop_var} vtkSurface vtkDiscreteRemeshing
	vtkCommon vtkFiltering vtkGraphics vtkImaging vtkIO vtkRendering vtkHybrid)
ENDFOREACH(loop_var)

## eof - CMakeLists.txt
