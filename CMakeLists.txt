cmake_minimum_required(VERSION 3.20)
set(NAME cfd-numerical-solver)
project(${NAME} LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

#set(CMAKE_CXX_COMPILER clang++)

# Подключаем модуль FetchContent
include(FetchContent)

# Загружаем xtl (базовая зависимость для xtensor)
FetchContent_Declare(
        xtl
        GIT_REPOSITORY https://github.com/xtensor-stack/xtl.git
        GIT_TAG 0.8.1
        GIT_SHALLOW TRUE
)

# Загружаем xtensor
FetchContent_Declare(
        xtensor
        GIT_REPOSITORY https://github.com/xtensor-stack/xtensor.git
        GIT_TAG 0.27.1
        GIT_SHALLOW TRUE
)

# Загружаем yaml-cpp
FetchContent_Declare(
        yaml-cpp
        GIT_REPOSITORY https://github.com/jbeder/yaml-cpp.git
        GIT_TAG 0.8.0
        GIT_SHALLOW TRUE
)

# Делаем зависимости доступными
FetchContent_MakeAvailable(xtl xtensor yaml-cpp)

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

target_link_libraries(${NAME} PRIVATE xtensor yaml-cpp::yaml-cpp)
