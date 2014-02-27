PROJECT(vtkSurface)

# --------------------- OPTIONS OF THE PROJECT -----------------------------
# Build static lib by default

OPTION(BUILD_SHARED_LIBS "Build vtkSurface with shared libraries." ON)
OPTION(BUILD_DOCUMENTATION "Build the html documentation with doxygen." OFF)

OPTION(BUILD_DISCRETEREMESHING "Build Discrete Remeshing Tools." ON)
OPTION(BUILD_VOLUMEPROCESSING "Build Volume tools" ON)

OPTION(USE_MULTITHREADING "Use Multithreaded classes" OFF)

OPTION(USE_LLOYDCLUSTERING "Use Lloyd relaxations instead of new algorithm for clustering" OFF)
MARK_AS_ADVANCED(USE_LLOYDCLUSTERING)

#
SET (EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin CACHE PATH "Single output directory for building all executables.")
SET (LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin CACHE PATH "Single output directory for building all libraries.")
MARK_AS_ADVANCED(LIBRARY_OUTPUT_PATH EXECUTABLE_OUTPUT_PATH)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/ACVD.json
    ${EXECUTABLE_OUTPUT_PATH}/ACVD.json COPYONLY)



# -------------------- DEPENDENCIES TO EXTERNALS ---------------------------
# Find and load VTK settings.


cmake_minimum_required(VERSION 2.8.7)
find_package(VTK REQUIRED NO_MODULE)
include(${VTK_USE_FILE})

SET(VTKSURFACE_PROJECT_DIR ${PROJECT_SOURCE_DIR})
SET(VTKSURFACE_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/Common)
SET(VTKDISCRETEREMESHING_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/DiscreteRemeshing)
SET(VTKVOLUMEPROCESSING_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/VolumeProcessing)
# --------------------------------------------------------------------------
# Find and load Doxygen settings (when required)


IF(BUILD_DOCUMENTATION)
    FIND_PACKAGE(Doxygen)
ENDIF(BUILD_DOCUMENTATION)

# --------------------------------------------------------------------
# Switch to enable Multithreaded clustering

IF(USE_MULTITHREADING)
    ADD_DEFINITIONS (-D DOmultithread)
ENDIF(USE_MULTITHREADING)
IF(USE_LLOYDCLUSTERING)
    ADD_DEFINITIONS (-D DOLloydClustering)
ENDIF(USE_LLOYDCLUSTERING)

# -------------------------- SUBDIRS TO HANDLE -----------------------------

# The main library
ADD_SUBDIRECTORY(Common)

IF(BUILD_DISCRETEREMESHING)
    ADD_SUBDIRECTORY(DiscreteRemeshing)
ENDIF(BUILD_DISCRETEREMESHING)

IF(BUILD_VOLUMEPROCESSING)
    ADD_SUBDIRECTORY(VolumeProcessing)
ENDIF(BUILD_VOLUMEPROCESSING)

IF(BUILD_DOCUMENTATION)
    ADD_SUBDIRECTORY(doc)
ENDIF(BUILD_DOCUMENTATION)

## eof - CMakeLists.txt
