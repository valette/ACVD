#-----------------------------------------------------------------------------
# Build the documentation
set(PROJECT_HTML_OUTPUT "html.user")
set(PROJECT_DOC_LOGFILE "${PROJECT_BINARY_DIR}/doc/${PROJECT_HTML_OUTPUT}.log")

# Doc input
set(PROJECT_DOC_INPUT 
   ${PROJECT_SOURCE_DIR}/doc/DoxyMainPage.txt
)

set(PROJECT_DOC_INPUT 
${PROJECT_SOURCE_DIR}/Common
${PROJECT_SOURCE_DIR}/DiscreteRemeshing
${PROJECT_SOURCE_DIR}/VolumeProcessing
)
string(REGEX REPLACE ";" " " PROJECT_DOC_INPUT "${PROJECT_DOC_INPUT}")

# Doc exclude
set(PROJECT_DOC_EXCLUDE "")
string(REGEX REPLACE ";" " " PROJECT_DOC_EXCLUDE "${PROJECT_DOC_EXCLUDE}")

#-----------------------------------------------------------------------------
# DOT verification
if(DOT)
   get_filename_component(PROJECT_DOC_DOT_PATH ${DOT} PATH)
   set(PROJECT_DOC_HAVE_DOT "YES")
else(DOT)
   set(PROJECT_DOC_DOT_PATH "")
   set(PROJECT_DOC_HAVE_DOT "NO")
endif(DOT)

#-----------------------------------------------------------------------------
# Create file and project
configure_file(
   ${PROJECT_SOURCE_DIR}/doc/doxygen.config.in
   ${PROJECT_BINARY_DIR}/doc/DoxyfileUsers
   IMMEDIATE
)

add_custom_target(doc-user
   ALL
   ${DOXYGEN}
   ${PROJECT_BINARY_DIR}/doc/DoxyfileUsers
)

install( DIRECTORIES ${PROJECT_HTML_OUTPUT}
  DESTINATION ${INSTALL_DOC_DIR}
)
#-----------------------------------------------------------------------------
