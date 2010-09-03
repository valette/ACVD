



#-----------------------------------------------------------------------------
# INSTALLATION OF FindLibrary.cmake in CMake/Modules dir

# ON WINDOWS : NO INSTALL MECHANISM; DIRTY COPY ALLOWED 
IF (WIN32)
 # CONFIGURE_FILE(
 #   ${INSTALL_LIBRARY_FILES_DIR}/FindLibrary.cmake.in
 #   ${CMAKE_ROOT}/Modules/Find${LIBRARY_NAME}.cmake
 #   @ONLY IMMEDIATE
 #   )
 # CONFIGURE_FILE(
 #   ${INSTALL_LIBRARY_FILES_DIR}/UseLibrary.cmake.in
 #   ${CURRENT_BINARY_DIR}/Use${LIBRARY_NAME}.cmake
 #   @ONLY IMMEDIATE
 #   )
 # # Create the LibraryConfig.cmake file containing the Library configuration.
 # INCLUDE (${INSTALL_LIBRARY_FILES_DIR}/GenerateLibraryConfig.cmake)
 # # Save the compiler settings so another project can import them.
 # INCLUDE(${CMAKE_ROOT}/Modules/CMakeExportBuildSettings.cmake)
 # CMAKE_EXPORT_BUILD_SETTINGS(${CMAKE_CURRENT_BINARY_DIR}/${LIBRARY_NAME}BuildSettings.cmake)

#  INSTALL(
#    FILES ${CMAKE_CURRENT_BINARY_DIR}/${LIBRARY_NAME}BuildSettings.cmake
#    DESTINATION lib/${LIBRARY_NAME}
#    )

  # ON LINUX : FILES MUST BE INSTALLED BY make install 
ELSE (WIN32)

  CONFIGURE_FILE(
    ${INSTALL_LIBRARY_FILES_DIR}/FindLibrary.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/Find${LIBRARY_NAME}.cmake
    @ONLY IMMEDIATE
    )
  INSTALL( 
    FILES ${CMAKE_CURRENT_BINARY_DIR}/Find${LIBRARY_NAME}.cmake
    DESTINATION ${CMAKE_ROOT}/Modules 
    )
  CONFIGURE_FILE(
    ${INSTALL_LIBRARY_FILES_DIR}/UseLibrary.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/Use${LIBRARY_NAME}.cmake
    @ONLY IMMEDIATE
    )
  INSTALL(
    FILES ${CMAKE_CURRENT_BINARY_DIR}/Use${LIBRARY_NAME}.cmake
    DESTINATION lib/${LIBRARY_NAME}
    )

  # Create the LibraryConfig.cmake file containing the lib configuration.
  INCLUDE (${INSTALL_LIBRARY_FILES_DIR}/GenerateLibraryConfig.cmake)

  # Save the compiler settings so another project can import them.
  INCLUDE(${CMAKE_ROOT}/Modules/CMakeExportBuildSettings.cmake)
  CMAKE_EXPORT_BUILD_SETTINGS(${CMAKE_CURRENT_BINARY_DIR}/${LIBRARY_NAME}BuildSettings.cmake)
  
  INSTALL(
    FILES ${CMAKE_CURRENT_BINARY_DIR}/${LIBRARY_NAME}BuildSettings.cmake
    DESTINATION lib/${LIBRARY_NAME}
    )

ENDIF(WIN32)
#-----------------------------------------------------------------------------
