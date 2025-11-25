include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/DoxygenFetcher.cmake)

# Handle Doxygen
setup_doxygen()

if(DOXYGEN_FOUND)
    # Fetch doxygen-awesome-css theme only if Doxygen is available
    FetchContent_Declare(
        doxygen-awesome-css
        URL https://github.com/jothepro/doxygen-awesome-css/archive/refs/heads/main.zip
    )
    FetchContent_MakeAvailable(doxygen-awesome-css)

    # Get the path where doxygen-awesome-css was downloaded
    FetchContent_GetProperties(doxygen-awesome-css SOURCE_DIR AWESOME_CSS_DIR)

    # Set paths for doxygen-awesome-css files
    set(DOXYGEN_AWESOME_CSS "${AWESOME_CSS_DIR}/doxygen-awesome.css")
    set(DOXYGEN_AWESOME_DARKMODE_JS "${AWESOME_CSS_DIR}/doxygen-awesome-darkmode-toggle.js")
    set(DOXYGEN_AWESOME_COPY_JS "${AWESOME_CSS_DIR}/doxygen-awesome-fragment-copy-button.js")
    set(DOXYGEN_AWESOME_PARAGRAPH_JS "${AWESOME_CSS_DIR}/doxygen-awesome-paragraph-link.js")

    # Configure Doxyfile
    set(DOXYFILE_IN ${CMAKE_CURRENT_SOURCE_DIR}/doc/Doxyfile.in)
    set(DOXYFILE_OUT ${CMAKE_CURRENT_BINARY_DIR}/doc/Doxyfile)
    configure_file(${DOXYFILE_IN} ${DOXYFILE_OUT} @ONLY)

    # Copy custom files to build directory
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/doc/custom_header.html
                   ${CMAKE_CURRENT_BINARY_DIR}/doc/custom_header.html COPYONLY)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/doc/custom_footer.html
                   ${CMAKE_CURRENT_BINARY_DIR}/doc/custom_footer.html COPYONLY)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/doc/extra_stylesheet.css
                   ${CMAKE_CURRENT_BINARY_DIR}/doc/extra_stylesheet.css COPYONLY)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/doc/DoxygenLayout.xml
                   ${CMAKE_CURRENT_BINARY_DIR}/doc/DoxygenLayout.xml COPYONLY)

    # Copy logo and favicon if they exist (deduplicated)
    if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/doc/full_logo.svg)
        configure_file(${CMAKE_CURRENT_SOURCE_DIR}/doc/full_logo.svg
                       ${CMAKE_CURRENT_BINARY_DIR}/doc/full_logo.svg COPYONLY)
    endif()
    if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/doc/favicon.ico)
        configure_file(${CMAKE_CURRENT_SOURCE_DIR}/doc/favicon.ico
                       ${CMAKE_CURRENT_BINARY_DIR}/doc/favicon.ico COPYONLY)
    endif()

    # Add custom target to generate documentation
    add_custom_target(docs
        COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYFILE_OUT}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating documentation with Doxygen"
        VERBATIM
    )

    message(STATUS "Doxygen configured - use 'make docs' or 'cmake --build . --target docs' to generate documentation")
endif()
