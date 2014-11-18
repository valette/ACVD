project(ACVD)

cmake_minimum_required(VERSION 2.8.7)
find_package(VTK REQUIRED NO_MODULE)
include(${VTK_USE_FILE})

option(BUILD_SHARED_LIBS "Build vtkSurface with shared libraries." ON)
option(BUILD_DOCUMENTATION "Build the html documentation with doxygen." OFF)

option(BUILD_DISCRETEREMESHING "Build Discrete Remeshing Tools." ON)
option(BUILD_VOLUMEPROCESSING "Build Volume tools" ON)

option(BUILD_EXAMPLES "Build examples" ON)

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin CACHE PATH "Single output directory for building all executables.")

set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin CACHE PATH "Single output directory for building all libraries.")
mark_as_advanced(LIBRARY_OUTPUT_PATH EXECUTABLE_OUTPUT_PATH)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/ACVD.json ${EXECUTABLE_OUTPUT_PATH}/ACVD.json COPYONLY)

set(VTKSURFACE_PROJECT_DIR ${PROJECT_SOURCE_DIR})
set(VTKSURFACE_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/Common)
set(VTKDISCRETEREMESHING_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/DiscreteRemeshing)
set(VTKVOLUMEPROCESSING_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/VolumeProcessing)

if(BUILD_DOCUMENTATION)
    find_package(Doxygen)
    add_subdirectory(doc)
endif(BUILD_DOCUMENTATION)

add_subdirectory(Common)

if(BUILD_DISCRETEREMESHING)
    add_subdirectory(DiscreteRemeshing)
endif(BUILD_DISCRETEREMESHING)

if(BUILD_VOLUMEPROCESSING)
    add_subdirectory(VolumeProcessing)
endif(BUILD_VOLUMEPROCESSING)

