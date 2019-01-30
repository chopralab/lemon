find_package(Sphinx REQUIRED)
find_package(Doxygen REQUIRED)

execute_process(
    COMMAND ${DOXYGEN_EXECUTABLE} --version
    OUTPUT_VARIABLE DOXYGEN_VERSION
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

if ("${DOXYGEN_VERSION}" VERSION_LESS 1.8.13)
    message(WARNING
        "Your version of doxygen is too old (${DOXYGEN_VERSION}).\n"
        "Please upgrade to doxygen >= 1.8.13 for optimal output."
    )
endif()

configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in
    ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
    @ONLY
)

configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/conf.py.in
    ${CMAKE_CURRENT_BINARY_DIR}/conf.py
    @ONLY
)

file(GLOB_RECURSE HEADERS ${PROJECT_SOURCE_DIR}/include/*.hpp)

add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/doxygen/xml/index.xml
    COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
    DEPENDS ${HEADERS} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
    COMMENT "Generating API documentation with Doxygen"
)

add_custom_target(doxygen ALL
    DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/doxygen/xml/index.xml
)

set(SPHINX_CONF_DIR ${CMAKE_CURRENT_BINARY_DIR})
set(SPHINX_INPUT_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)
set(SPHINX_DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
set(SPHINX_DEPENDENCIES doxygen)

set(sphinx_target_base_name doc)
Sphinx_add_targets(${sphinx_target_base_name} ${SPHINX_CONF_DIR} ${SPHINX_INPUT_DIR} ${SPHINX_DESTINATION} ${SPHINX_DEPENDENCIES})

if( ${SPHINX_HTML_OUTPUT} )
    install( DIRECTORY ${SPHINX_DESTINATION}/html
        DESTINATION share/lemon/doc
        COMPONENT Doc
        PATTERN "${SPHINX_DESTINATION}/html/*"
    )
endif()
