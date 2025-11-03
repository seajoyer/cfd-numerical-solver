cmake_minimum_required(VERSION 3.20)
set(NAME cfd-numerical-solver)
project(${NAME} LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(xtensor REQUIRED)

include_directories(${PROJECT_SOURCE_DIR}/include)

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/config/*.cpp
        ${PROJECT_SOURCE_DIR}/src/data/*.cpp
        ${PROJECT_SOURCE_DIR}/src/solver/*.cpp
        ${PROJECT_SOURCE_DIR}/src/visualization/*.cpp
        ${PROJECT_SOURCE_DIR}/src/bc/*.cpp
        ${PROJECT_SOURCE_DIR}/src/*.cpp)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/config/*.hpp
        ${PROJECT_SOURCE_DIR}/include/data/*.hpp
        ${PROJECT_SOURCE_DIR}/include/solver/*.hpp
        ${PROJECT_SOURCE_DIR}/include/visualization/*.hpp
        ${PROJECT_SOURCE_DIR}/include/bc/*.hpp
        ${PROJECT_SOURCE_DIR}/include/*.hpp)

add_executable(${NAME} src/main.cpp ${sources} ${headers})

target_link_libraries(${NAME} PRIVATE xtensor)
