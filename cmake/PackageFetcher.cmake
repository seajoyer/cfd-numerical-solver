macro(fetch_package target_name repo tag)
    if(NOT TARGET ${target_name} AND NOT TARGET ${target_name}::${target_name})
        find_package(${target_name} QUIET CONFIG)
        if(NOT ${target_name}_FOUND AND NOT TARGET ${target_name}::${target_name})
            message(STATUS "[fetch_package] ${target_name} not found → downloading ${tag}")

            FetchContent_Declare(
                ${target_name}
                GIT_REPOSITORY ${repo}
                GIT_TAG        ${tag}
                GIT_SHALLOW    TRUE
                FIND_PACKAGE_ARGS CONFIG
            )

            # special handling for VTK
            if(${target_name} STREQUAL "VTK")
                message(STATUS "[fetch_package] Configuring minimal VTK build (CommonCore + CommonDataModel + IOCore + IOLegacy only)...")

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
            endif()

            FetchContent_MakeAvailable(${target_name})
        else()
            message(STATUS "[fetch_package] Found ${target_name}")
        endif()
    endif()
endmacro()
