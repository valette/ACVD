# Generate the ${LIBRARY_NAME}Config.cmake file in the build tree. 
# Also configure one for installation. 
# The file tells external projects how to use ${LIBRARY_NAME}.

#-----------------------------------------------------------------------------
# Global settings
SET(LIB_MAJOR_VERSION_CONFIG ${${LIBRARY_NAME}_MAJOR_VERSION})
SET(LIB_MINOR_VERSION_CONFIG ${${LIBRARY_NAME}_MINOR_VERSION})
SET(LIB_BUILD_VERSION_CONFIG ${${LIBRARY_NAME}_BUILD_VERSION})
SET(LIB_VERSION_CONFIG 
  ${${LIBRARY_NAME}_MAJOR_VERSION}.${${LIBRARY_NAME}_MINOR_VERSION}.${${LIBRARY_NAME}_BUILD_VERSION})

SET(LIB_REQUIRED_C_FLAGS_CONFIG ${${LIBRARY_NAME}_REQUIRED_C_FLAGS})
SET(LIB_REQUIRED_CXX_FLAGS_CONFIG ${${LIBRARY_NAME}_CXX_FLAGS})
SET(LIB_REQUIRED_LINK_FLAGS_CONFIG ${${LIBRARY_NAME}_LINK_FLAGS})
#-----------------------------------------------------------------------------


#-----------------------------------------------------------------------------
# Settings specific to the install tree.

# The "use" file.
SET(LIB_USE_FILE_CONFIG 
  ${CMAKE_INSTALL_PREFIX}/lib/${LIBRARY_NAME}/Use${LIBRARY_NAME}.cmake)

# The build settings file.
SET(LIB_BUILD_SETTINGS_FILE_CONFIG
  ${CMAKE_INSTALL_PREFIX}/lib/${LIBRARY_NAME}/${LIBRARY_NAME}BuildSettings.cmake)

# Include directories.
SET(LIB_INCLUDE_DIRS_CONFIG
  ${CMAKE_INSTALL_PREFIX}/include/${LIBRARY_NAME}
  )

# Link directories.
SET(LIB_LIBRARY_DIRS_CONFIG ${CMAKE_INSTALL_PREFIX}/lib)

#-----------------------------------------------------------------------------
# Configure ${LIBRARY_NAME}Config.cmake for the install tree.
CONFIGURE_FILE(${INSTALL_LIBRARY_FILES_DIR}/LibraryConfig.cmake.in
               ${PROJECT_BINARY_DIR}/cmake/${LIBRARY_NAME}Config.cmake 
	       @ONLY IMMEDIATE)

INSTALL(
  FILES ${PROJECT_BINARY_DIR}/cmake/${LIBRARY_NAME}Config.cmake
  DESTINATION lib/${LIBRARY_NAME}
  )
#INSTALL_FILES(/lib/${LIBRARY_NAME} .cmake ${LIBRARY_NAME}Config)

#-----------------------------------------------------------------------------
# Settings specific to the build tree.

# The "use" file.
SET(LIB_USE_FILE_CONFIG 
  ${CMAKE_CURRENT_BINARY_DIR}/Use${LIBRARY_NAME}.cmake)

# The build settings file.
SET(LIB_BUILD_SETTINGS_FILE_CONFIG 
  ${CMAKE_CURRENT_BINARY_DIR}/${LIBRARY_NAME}BuildSettings.cmake)

# Library directory.
SET(LIB_LIBRARY_DIRS_CONFIG ${LIBRARY_OUTPUT_PATH})

# Determine the include directories needed.
SET(LIB_INCLUDE_DIRS_CONFIG
  ${CMAKE_CURRENT_SOURCE_DIR}
)

#-----------------------------------------------------------------------------
# Configure ${LIBRARY_NAME}Config.cmake for the build tree.
CONFIGURE_FILE(${INSTALL_LIBRARY_FILES_DIR}/LibraryConfig.cmake.in
               ${LIBRARY_NAME_BUILD_TREE_CONFIG}/${LIBRARY_NAME}Config.cmake 
	       @ONLY IMMEDIATE)
