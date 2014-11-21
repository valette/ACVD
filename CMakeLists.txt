project(ACVD)

cmake_minimum_required(VERSION 2.8.7)

if( POLICY CMP0042 )
  cmake_policy( SET CMP0042 NEW )
endif()

# version
set( ACVD_MAJOR_VERSION 1 )
set( ACVD_MINOR_VERSION 0 )
set( ACVD_PATCH_VERSION 0 )
set( ACVD_VERSION ${ACVD_MAJOR_VERSION}.${ACVD_MINOR_VERSION}.${ACVD_PATCH_VERSION} )

set( ACVD_LIB_VERSION ${ACVD_VERSION} )
set( ACVD_LIB_SOVERSION ${ACVD_MAJOR_VERSION}.${ACVD_MINOR_VERSION} )

# installation
set(INSTALL_LIB_DIR lib CACHE PATH "Installation directory for libraries")
set(INSTALL_BIN_DIR bin CACHE PATH "Installation directory for executables")
set(INSTALL_INCLUDE_DIR include CACHE PATH "Installation directory for header files")
set(INSTALL_DOC_DIR share/doc/ACVD CACHE PATH "Installation directory for html doc" )
set( INSTALL_CMAKE_DIR lib/cmake/ACVD "Installation directory for CMake files" )

# find_package
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

