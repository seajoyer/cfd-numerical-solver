cmake_minimum_required(VERSION 3.25)

set(NAME cfd-numerical-solver)
project(${NAME} LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

find_package(xtensor REQUIRED)
find_package(VTK REQUIRED)

include(FetchContent)

# Set to prefer system packages
set(FETCHCONTENT_TRY_FIND_PACKAGE_MODE ALWAYS)

macro(safe_fetch_content package_name repo tag)
    if(NOT TARGET ${package_name})
        find_package(${package_name} QUIET)
        if(NOT ${package_name}_FOUND)
            FetchContent_Declare(
                ${package_name}
                GIT_REPOSITORY ${repo}
                GIT_TAG ${tag}
                GIT_SHALLOW TRUE
                FIND_PACKAGE_ARGS NAMES ${package_name}
            )
            FetchContent_MakeAvailable(${package_name})
        endif()
    endif()
endmacro()

# Use the macro for each dependency
safe_fetch_content(xtl https://github.com/xtensor-stack/xtl.git 0.8.1)
safe_fetch_content(xtensor https://github.com/xtensor-stack/xtensor.git 0.27.1)
safe_fetch_content(yaml-cpp https://github.com/jbeder/yaml-cpp.git 0.8.0)
safe_fetch_content(cxxopts https://github.com/jarro2783/cxxopts.git 3.3.1)

find_package(Doxygen REQUIRED)

if(DOXYGEN_FOUND)
    # Fetch doxygen-awesome-css theme
    include(FetchContent)
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

    # Copy logo if it exists
    if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/doc/full_logo.svg)
        configure_file(${CMAKE_CURRENT_SOURCE_DIR}/doc/full_logo.svg
                       ${CMAKE_CURRENT_BINARY_DIR}/doc/full_logo.svg COPYONLY)
    endif()
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

# Explicitly list your sources
add_executable(${NAME}
    src/main.cpp

    src/Simulation.cpp

    src/solver/SolverFactory.cpp
    src/solver/GodunovSolver.cpp
    src/solver/GodunovKolganRodionovSolver.cpp
    src/solver/AnalyticalSolver.cpp
    include/solver/Solver.hpp
    src/solver/EOS.cpp
    src/solver/TimeStepCalculator.cpp
    src/solver/PositivityLimiter.cpp

    src/reconstruction/P0Reconstruction.cpp
    src/reconstruction/P1Reconstruction.cpp
    include/reconstruction/Reconstruction.hpp

    src/riemann/HLLRiemannSolver.cpp
    src/riemann/HLLCRiemannSolver.cpp
    src/riemann/ExactIdealGasRiemannSolver.cpp
    src/riemann/AcousticRiemannSolver.cpp

    src/output/WriterFactory.cpp
    src/output/VTKWriter.cpp

    src/data/DataLayer.cpp
    src/data/Variables.cpp

    src/config/ConfigParser.cpp
    include/config/Settings.hpp
    include/config/InitialConditions.hpp

    src/bc/BoundaryFactory.cpp
    src/bc/BoundaryManager.cpp
    src/bc/FreeStreamBoundary.cpp
    src/bc/InletBoundary.cpp
    src/bc/OutletBoundary.cpp
    src/bc/PeriodicBoundary.cpp
    src/bc/ReflectiveBoundary.cpp
    src/bc/NonReflectiveBoundary.cpp
    src/bc/SymmetryBoundary.cpp
    src/bc/WallBoundary.cpp
    include/bc/BoundaryCondition.hpp

    src/utils/StringUtils.cpp
)

target_include_directories(${NAME}
    PRIVATE
        ${PROJECT_SOURCE_DIR}/include
)

target_link_libraries(${NAME}
    PRIVATE
        xtensor
        yaml-cpp::yaml-cpp
        cxxopts::cxxopts
        ${VTK_LIBRARIES}
)
