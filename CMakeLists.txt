project(ACVD)

cmake_minimum_required(VERSION 3.0.0)

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
set(INSTALL_CMAKE_DIR lib/cmake/ACVD CACHE PATH "Installation directory for CMake files" )

# Make relative paths absolute (needed later on)
foreach(p LIB BIN INCLUDE CMAKE)
  set(var INSTALL_${p}_DIR)
  if(NOT IS_ABSOLUTE "${${var}}")
    set(${var} "${CMAKE_INSTALL_PREFIX}/${${var}}")
  endif()
endforeach()

# find_package
find_package(VTK REQUIRED NO_MODULE)
if(VTK_USE_FILE)
  include(${VTK_USE_FILE})
else()
  include_directories(${VTK_INCLUDE_DIRS})
endif()

option(BUILD_SHARED_LIBS "Build vtkSurface with shared libraries." ON)
option(BUILD_DOCUMENTATION "Build the html documentation with doxygen." OFF)

option(BUILD_DISCRETEREMESHING "Build Discrete Remeshing Tools." ON)
option(BUILD_VOLUMEPROCESSING "Build Volume tools" ON)

option(BUILD_EXAMPLES "Build examples" ON)

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin CACHE PATH "Single output directory for building all executables.")

set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin CACHE PATH "Single output directory for building all libraries.")
mark_as_advanced(LIBRARY_OUTPUT_PATH EXECUTABLE_OUTPUT_PATH)

#configure_file(${CMAKE_CURRENT_SOURCE_DIR}/ACVD.json ${EXECUTABLE_OUTPUT_PATH}/ACVD.json COPYONLY)

set(VTKSURFACE_PROJECT_DIR ${PROJECT_SOURCE_DIR})
set(VTKSURFACE_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/Common)
set(VTKDISCRETEREMESHING_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/DiscreteRemeshing)
set(VTKVOLUMEPROCESSING_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/VolumeProcessing)

if(BUILD_DOCUMENTATION)
    find_package(Doxygen)
    add_subdirectory(doc)
endif(BUILD_DOCUMENTATION)

add_subdirectory(Common)
set( ACVD_LIBRARIES vtkSurface )

if(BUILD_DISCRETEREMESHING)
    add_subdirectory(DiscreteRemeshing)
    set( ACVD_LIBRARIES ${ACVD_LIBRARIES} vtkDiscreteRemeshing )
endif(BUILD_DISCRETEREMESHING)

if(BUILD_VOLUMEPROCESSING)
    add_subdirectory(VolumeProcessing)
    set( ACVD_LIBRARIES ${ACVD_LIBRARIES} vtkVolumeProcessing )
endif(BUILD_VOLUMEPROCESSING)

export( TARGETS ${_targets}
  FILE "${ACVD_BINARY_DIR}/ACVDTargets.cmake"
)

# Export the package for use from the build-tree
# (this registers the build-tree with a global CMake-registry)
export(PACKAGE ACVD)

# -------------------------------------------------------------
# for the build tree
set( CONF_INCLUDE_DIRS
  "${ACVD_SOURCE_DIR}/Common"
)

if( BUILD_DISCRETEREMESHING )
  set( CONF_INCLUDE_DIRS "${CONF_INCLUDE_DIRS}" "${ACVD_SOURCE_DIR}/DiscreteRemeshing" )
endif()

if( BUILD_VOLUMEPROCESSING )
  set( CONF_INCLUDE_DIRS "${CONF_INCLUDE_DIRS}" "${ACVD_SOURCE_DIR}/VolumeProcessing" )
endif()

configure_file( ${ACVD_SOURCE_DIR}/ACVDConfig.cmake.in
  ${ACVD_BINARY_DIR}/ACVDConfig.cmake
  @ONLY
)
# for the install tree
file(RELATIVE_PATH REL_INCLUDE_DIR "${INSTALL_CMAKE_DIR}"
   "${INSTALL_INCLUDE_DIR}")

set( CONF_INCLUDE_DIRS
  "\${ACVD_CMAKE_DIR}/${REL_INCLUDE_DIR}/ACVD/Common"
)
if( BUILD_DISCRETEREMESHING )
  set( CONF_INCLUDE_DIRS "${CONF_INCLUDE_DIRS}" "\${ACVD_CMAKE_DIR}/${REL_INCLUDE_DIR}/ACVD/DiscreteRemeshing" )
endif()

if( BUILD_VOLUMEPROCESSING )
  set( CONF_INCLUDE_DIRS "${CONF_INCLUDE_DIRS}" "\${ACVD_CMAKE_DIR}/${REL_INCLUDE_DIR}/ACVD/VolumeProcessing" )
endif()

configure_file( ${ACVD_SOURCE_DIR}/ACVDConfig.cmake.in
  ${ACVD_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/ACVDConfig.cmake
  @ONLY
)
configure_file( ${ACVD_SOURCE_DIR}/ACVDConfigVersion.cmake.in
  ${ACVD_BINARY_DIR}/ACVDConfigVersion.cmake
  @ONLY
)
install( FILES
    "${ACVD_BINARY_DIR}/ACVDConfig.cmake"
    "${ACVD_BINARY_DIR}/ACVDConfigVersion.cmake"
  DESTINATION ${INSTALL_CMAKE_DIR} COMPONENT dev
)
install(EXPORT ACVDTargets DESTINATION "${INSTALL_CMAKE_DIR}" COMPONENT dev)

