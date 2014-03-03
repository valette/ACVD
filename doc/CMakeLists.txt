#-----------------------------------------------------------------------------
# Build the documentation
SET(PROJECT_HTML_OUTPUT "html.user")
SET(PROJECT_DOC_LOGFILE "${PROJECT_BINARY_DIR}/doc/${PROJECT_HTML_OUTPUT}.log")

# Doc input
SET(PROJECT_DOC_INPUT 
   ${PROJECT_SOURCE_DIR}/doc/DoxyMainPage.txt
)

SET(PROJECT_DOC_INPUT 
${PROJECT_SOURCE_DIR}/Common
${PROJECT_SOURCE_DIR}/DiscreteRemeshing
${PROJECT_SOURCE_DIR}/VolumeProcessing
)
STRING(REGEX REPLACE ";" " " PROJECT_DOC_INPUT "${PROJECT_DOC_INPUT}")

# Doc exclude
SET(PROJECT_DOC_EXCLUDE "")
STRING(REGEX REPLACE ";" " " PROJECT_DOC_EXCLUDE "${PROJECT_DOC_EXCLUDE}")

#-----------------------------------------------------------------------------
# DOT verification
IF(DOT)
   GET_FILENAME_COMPONENT(PROJECT_DOC_DOT_PATH ${DOT} PATH)
   SET(PROJECT_DOC_HAVE_DOT "YES")
ELSE(DOT)
   SET(PROJECT_DOC_DOT_PATH "")
   SET(PROJECT_DOC_HAVE_DOT "NO")
ENDIF(DOT)

#-----------------------------------------------------------------------------
# Create file and project
CONFIGURE_FILE(
   ${PROJECT_SOURCE_DIR}/doc/doxygen.config.in
   ${PROJECT_BINARY_DIR}/doc/DoxyfileUsers
   IMMEDIATE
)

ADD_CUSTOM_TARGET(doc-user
   ALL
   ${DOXYGEN}
   ${PROJECT_BINARY_DIR}/doc/DoxyfileUsers
)

#-----------------------------------------------------------------------------
