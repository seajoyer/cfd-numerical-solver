macro(fetch_package target_name repo tag)
    if (NOT TARGET ${target_name} AND NOT TARGET ${target_name}::${target_name})
        find_package(${target_name} QUIET CONFIG)
        if (NOT ${target_name}_FOUND AND NOT TARGET ${target_name}::${target_name})
            message(STATUS "[fetch_package] ${target_name} not found → downloading ${tag}")

            FetchContent_Declare(
                    ${target_name}
                    GIT_REPOSITORY ${repo}
                    GIT_TAG ${tag}
                    GIT_SHALLOW TRUE
                    FIND_PACKAGE_ARGS CONFIG
            )

            # special handling for VTK
            if (${target_name} STREQUAL "VTK")
                message(STATUS "[fetch_package] Configuring minimal VTK build (CommonCore + CommonDataModel + IOLegacy only)...")

                # Disable extras
                set(VTK_BUILD_TESTING OFF CACHE BOOL "")
                set(VTK_WRAP_PYTHON OFF CACHE BOOL "")
                set(VTK_WRAP_JAVA OFF CACHE BOOL "")
                set(VTK_USE_MPI OFF CACHE BOOL "")
                set(VTK_USE_TK OFF CACHE BOOL "")
                set(VTK_BUILD_DOCUMENTATION OFF CACHE BOOL "")

                # Disable all groups to DONT_WANT (prevents unwanted modules)
                set(VTK_GROUP_ENABLE_Imaging DONT_WANT CACHE STRING "")
                set(VTK_GROUP_ENABLE_MPI DONT_WANT CACHE STRING "")
                set(VTK_GROUP_ENABLE_Qt DONT_WANT CACHE STRING "")
                set(VTK_GROUP_ENABLE_Rendering DONT_WANT CACHE STRING "")
                set(VTK_GROUP_ENABLE_StandAlone DONT_WANT CACHE STRING "")
                set(VTK_GROUP_ENABLE_Views DONT_WANT CACHE STRING "")
                set(VTK_GROUP_ENABLE_Web DONT_WANT CACHE STRING "")

                # Explicitly enable only required modules
                set(VTK_MODULE_ENABLE_VTK_CommonCore YES CACHE STRING "")
                set(VTK_MODULE_ENABLE_VTK_CommonDataModel YES CACHE STRING "")
                set(VTK_MODULE_ENABLE_VTK_IOLegacy YES CACHE STRING "")
            endif ()

            FetchContent_MakeAvailable(${target_name})
        else ()
            message(STATUS "[fetch_package] Found ${target_name}")
        endif ()
    endif ()
endmacro()

include(FetchContent)

function(fetch_gmsh)
    if(TARGET gmsh::gmsh)
        message(STATUS "[fetch_gmsh] Gmsh target already exists: gmsh::gmsh")
        return()
    endif()

    find_package(gmsh QUIET CONFIG)
    if(gmsh_FOUND AND TARGET gmsh::gmsh)
        message(STATUS "[fetch_gmsh] Found system Gmsh: gmsh::gmsh")
        return()
    endif()

    message(STATUS "[fetch_gmsh] Gmsh not found -> downloading official SDK 4.15.1")

    FetchContent_Declare(
            gmsh_sdk
            URL https://gmsh.info/bin/Linux/gmsh-4.15.1-Linux64-sdk.tgz
            DOWNLOAD_EXTRACT_TIMESTAMP TRUE
    )

    FetchContent_GetProperties(gmsh_sdk)
    if(NOT gmsh_sdk_POPULATED)
        FetchContent_Populate(gmsh_sdk)
    endif()

    set(GMSH_SDK_ROOT "${gmsh_sdk_SOURCE_DIR}")

    set(GMSH_INCLUDE_DIR "${GMSH_SDK_ROOT}/include")
    set(GMSH_LIBRARY_DIR "${GMSH_SDK_ROOT}/lib")
    set(GMSH_LIBRARY     "${GMSH_LIBRARY_DIR}/libgmsh.so")

    if(NOT EXISTS "${GMSH_INCLUDE_DIR}/gmsh.h")
        message(FATAL_ERROR "[fetch_gmsh] gmsh.h not found in SDK include directory: ${GMSH_INCLUDE_DIR}")
    endif()

    if(NOT EXISTS "${GMSH_LIBRARY}")
        message(FATAL_ERROR "[fetch_gmsh] libgmsh.so not found in SDK lib directory: ${GMSH_LIBRARY_DIR}")
    endif()

    add_library(gmsh::gmsh SHARED IMPORTED GLOBAL)
    set_target_properties(gmsh::gmsh PROPERTIES
            IMPORTED_LOCATION "${GMSH_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${GMSH_INCLUDE_DIR}"
    )

    message(STATUS "[fetch_gmsh] Using SDK imported target gmsh::gmsh")
endfunction()