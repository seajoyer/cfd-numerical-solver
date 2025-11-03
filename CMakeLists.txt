cmake_minimum_required(VERSION 3.20)
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

include_directories(${PROJECT_SOURCE_DIR}/include)

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/config/*.cpp
        ${PROJECT_SOURCE_DIR}/src/data/*.cpp
        ${PROJECT_SOURCE_DIR}/src/solver/*.cpp
        ${PROJECT_SOURCE_DIR}/src/visualization/*.cpp)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/config/*.hpp
        ${PROJECT_SOURCE_DIR}/include/data/*.hpp
        ${PROJECT_SOURCE_DIR}/include/solver/*.hpp
        ${PROJECT_SOURCE_DIR}/include/visualization/*.hpp)

add_executable(${NAME} src/main.cpp ${sources} ${headers})

target_link_libraries(${NAME} PRIVATE xtensor yaml-cpp::yaml-cpp ${VTK_LIBRARIES})
