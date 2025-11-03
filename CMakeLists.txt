cmake_minimum_required(VERSION 4.0)
set(NAME cfd-numerical-solver)
project(${NAME} LANGUAGES CXX)

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

# Explicitly list your sources
add_executable(${NAME}
    src/main.cpp

    src/config/ConfigParser.cpp
    src/config/Config.cpp
    src/config/StepWriter.cpp

    src/data/DataLayer.cpp

    src/solver/Solver.cpp
    src/solver/GodunovSolver.cpp

    src/visualization/VTKWriter.cpp

    src/bc/BoundaryManager.cpp
    src/bc/FreeStreamBoundary.cpp
    src/bc/InletBoundary.cpp
    src/bc/OutletBoundary.cpp
    src/bc/PeriodicBoundary.cpp
    src/bc/ReflectiveBoundary.cpp
    src/bc/SymmetryBoundary.cpp
    src/bc/WallBoundary.cpp 
)

target_include_directories(${NAME}
    PRIVATE
        ${PROJECT_SOURCE_DIR}/include
)

target_link_libraries(${NAME}
    PRIVATE
        xtensor
        yaml-cpp::yaml-cpp
        ${VTK_LIBRARIES}
)
