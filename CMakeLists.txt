cmake_minimum_required(VERSION 3.0.0)
project(ACVD VERSION 1.0.1)

set(ACVD_LIB_SOVERSION ${ACVD_MAJOR_VERSION}.${ACVD_MINOR_VERSION} )

# Build Env Options
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)

# Installation Env Options
set(INSTALL_LIB_DIR lib CACHE PATH "Installation directory for libraries")
set(INSTALL_BIN_DIR bin CACHE PATH "Installation directory for executables")
set(INSTALL_INCLUDE_DIR include CACHE PATH "Installation directory for header files")
set(INSTALL_DOC_DIR share/doc/ACVD CACHE PATH "Installation directory for html doc" )
set(INSTALL_CMAKE_DIR lib/cmake/ACVD CACHE PATH "Installation directory for CMake files")

# Build Options
option(BUILD_DOCUMENTATION "Build the html documentation with doxygen." OFF)
option(BUILD_DISCRETEREMESHING "Build Discrete Remeshing Tools." ON)
option(BUILD_VOLUMEPROCESSING "Build Volume tools" ON)
option(BUILD_EXAMPLES "Build examples" ON)

# App manifest
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/ACVD.json ${EXECUTABLE_OUTPUT_PATH}/ACVD.json COPYONLY)

# find dependencies
find_package(VTK REQUIRED NO_MODULE)
include(${VTK_USE_FILE})

if(BUILD_DOCUMENTATION)
    find_package(Doxygen)
    add_subdirectory(doc)
endif(BUILD_DOCUMENTATION)

# Build vtkSurface
add_subdirectory(vtkSurface)

# Build vtkDiscreteRemeshing
if(BUILD_DISCRETEREMESHING)
    add_subdirectory(vtkDiscreteRemeshing)
endif(BUILD_DISCRETEREMESHING)

# Build vtkVolumeProcessing
if(BUILD_VOLUMEPROCESSING)
    add_subdirectory(vtkVolumeProcessing)
endif(BUILD_VOLUMEPROCESSING)


#### Package Installation ####
# Configure package configs
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    ACVDConfigVersion.cmake COMPATIBILITY SameMajorVersion
)

configure_package_config_file(
    ACVDConfig.cmake.in
    ACVDConfig.cmake
    INSTALL_DESTINATION ${INSTALL_CMAKE_DIR}
)

# Install package configs
install(
    FILES
        ${CMAKE_BINARY_DIR}/ACVDConfig.cmake
        ${CMAKE_BINARY_DIR}/ACVDConfigVersion.cmake
    DESTINATION ${INSTALL_CMAKE_DIR}
)

# Install exported targets file
install(
    EXPORT ACVDTargets
    NAMESPACE "ACVD::"
    DESTINATION ${INSTALL_CMAKE_DIR}
)

