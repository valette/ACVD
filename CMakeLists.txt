PROJECT(VTKSURFACE)
CMAKE_MINIMUM_REQUIRED(VERSION 2.4)
#MARK_AS_ADVANCED( FORCE CMAKE_BACKWARDS_COMPATIBILITY )

# for CMake 2.6 corrected behaviour (see "cmake --help-policy CMP0003")
IF(COMMAND cmake_policy AND ${CMAKE_MAJOR_VERSION} EQUAL 2 AND ${CMAKE_MINOR_VERSION} GREATER 4)
  CMAKE_POLICY(SET CMP0003 NEW)
  CMAKE_POLICY(SET CMP0005 NEW)
ENDIF(COMMAND cmake_policy AND ${CMAKE_MAJOR_VERSION} EQUAL 2 AND ${CMAKE_MINOR_VERSION} GREATER 4)
# --------------------------------------------------------------------------


# --------------------- OPTIONS OF THE PROJECT -----------------------------
# Build static lib by default

OPTION(BUILD_SHARED_LIBS "Build vtkSurface with shared libraries." ON)
OPTION(BUILD_DOCUMENTATION "Build the html documentation with doxygen." OFF)
OPTION(BUILD_DISCRETEREMESHING "Build Discrete Remeshing Tools." ON)
OPTION(USE_MULTITHREADING "Use Multithreaded classes" OFF)
OPTION(USE_LLOYDCLUSTERING "Use Lloyd relaxations instead of new algorithm for clustering" OFF)
MARK_AS_ADVANCED(USE_LLOYDCLUSTERING)

# -------------------- DEPENDENCIES TO EXTERNALS ---------------------------
# Find and load VTK settings.

#INCLUDE(${CMAKE_ROOT}/Modules/FindVTK.cmake)
FIND_PACKAGE(VTK REQUIRED)
IF(VTK_FOUND)
    INCLUDE(${VTK_USE_FILE})
ENDIF(VTK_FOUND)


SET(VTKSURFACE_PROJECT_DIR ${PROJECT_SOURCE_DIR})
SET(VTKSURFACE_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/Common)
SET(VTKDISCRETEREMESHING_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/DiscreteRemeshing)

SET (EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin CACHE PATH "Single output directory for building all executables.")
SET (LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin CACHE PATH "Single output directory for building all libraries.")
MARK_AS_ADVANCED(LIBRARY_OUTPUT_PATH EXECUTABLE_OUTPUT_PATH)

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
SUBDIRS(Common)

IF(BUILD_DISCRETEREMESHING)
    SUBDIRS(DiscreteRemeshing)
ENDIF(BUILD_DISCRETEREMESHING)

