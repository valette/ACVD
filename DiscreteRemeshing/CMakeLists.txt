INCLUDE_DIRECTORIES(
    ${VTKSURFACE_INCLUDE_DIR}
)

OPTION(USE_MULTITHREADING "Use Multithreaded classes" OFF)
MARK_AS_ADVANCED(USE_MULTITHREADING)

OPTION(USE_LLOYDCLUSTERING "Use Lloyd relaxations instead of new algorithm for clustering" OFF)
MARK_AS_ADVANCED(USE_LLOYDCLUSTERING)


# --------------------------------------------------------------------
# Switch to enable Multithreaded clustering

IF(USE_MULTITHREADING)
    ADD_DEFINITIONS (-D DOmultithread)
ENDIF(USE_MULTITHREADING)
IF(USE_LLOYDCLUSTERING)
    ADD_DEFINITIONS (-D DOLloydClustering)
ENDIF(USE_LLOYDCLUSTERING)
 
 
 # --------------------------------------------------------------------------
 # Library compilation
 ADD_LIBRARY(
vtkDiscreteRemeshing 
vtkTestUniformClustering.cxx
xmlParser.cpp
vtkManifoldSimplification.cxx
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

INCLUDE_DIRECTORIES(
    ${VTK_INCLUDE_DIR}
    ${VTKSURFACE_INCLUDE_DIR}
    ${VTKDISCRETEREMESHING_INCLUDE_DIR}
)

SET(DISCRETEREMESHING_EXAMPLES
ACVD
ACVDQ
AnisotropicRemeshing
AnisotropicRemeshingQ
VolumeAnalysis
ManifoldSimplification
)

FOREACH(loop_var ${DISCRETEREMESHING_EXAMPLES})
	ADD_EXECUTABLE(${loop_var} ${loop_var}.cxx)
	TARGET_LINK_LIBRARIES(${loop_var} vtkSurface vtkDiscreteRemeshing
	vtkCommon vtkFiltering vtkGraphics vtkImaging vtkIO vtkRendering vtkHybrid)
ENDFOREACH(loop_var)

## eof - CMakeLists.txt
