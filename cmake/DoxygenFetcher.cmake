macro(setup_doxygen)
    find_package(Doxygen)
    if(NOT DOXYGEN_FOUND)
        set(DOXYGEN_VERSION "1.15.0")
        set(DOXYGEN_TAG "Release_1_15_0")

        message(STATUS "[setup_doxygen] Doxygen not found → attempting pre-built binary download for v${DOXYGEN_VERSION}")

        # Platform-specific binary settings
        if(WIN32 AND CMAKE_SIZEOF_VOID_P EQUAL 8)
            set(DOXYGEN_URL "https://github.com/doxygen/doxygen/releases/download/${DOXYGEN_TAG}/doxygen-${DOXYGEN_VERSION}.windows.x64.bin.zip")
            set(DOXYGEN_ARCHIVE_EXT ".windows.x64.bin.zip")
            set(DOXYGEN_BINARY_NAME "doxygen.exe")
        elseif(UNIX AND NOT APPLE AND CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64")
            set(DOXYGEN_URL "https://github.com/doxygen/doxygen/releases/download/${DOXYGEN_TAG}/doxygen-${DOXYGEN_VERSION}.linux.bin.tar.gz")
            set(DOXYGEN_ARCHIVE_EXT ".linux.bin.tar.gz")
            set(DOXYGEN_BINARY_NAME "doxygen")
        elseif(APPLE)
            set(DOXYGEN_URL "https://github.com/doxygen/doxygen/releases/download/${DOXYGEN_TAG}/Doxygen-${DOXYGEN_VERSION}.dmg")
            set(DOXYGEN_ARCHIVE_EXT ".dmg")
            set(DOXYGEN_BINARY_NAME "doxygen")
        else()
            message(WARNING "Pre-built Doxygen binary not available for this platform (${CMAKE_SYSTEM_NAME} ${CMAKE_SYSTEM_PROCESSOR}). Documentation target will not be available.")
            return()
        endif()

        # Download the binary archive
        set(DOXYGEN_DOWNLOAD_DIR "${CMAKE_BINARY_DIR}/doxygen-download")
        file(MAKE_DIRECTORY ${DOXYGEN_DOWNLOAD_DIR})
        set(DOXYGEN_ARCHIVE_FILE "${DOXYGEN_DOWNLOAD_DIR}/doxygen-${DOXYGEN_VERSION}${DOXYGEN_ARCHIVE_EXT}")

        file(DOWNLOAD ${DOXYGEN_URL} ${DOXYGEN_ARCHIVE_FILE}
             SHOW_PROGRESS
             STATUS DOXYGEN_DOWNLOAD_STATUS
        )
        list(GET DOXYGEN_DOWNLOAD_STATUS 0 DOXYGEN_DOWNLOAD_CODE)
        if(NOT DOXYGEN_DOWNLOAD_CODE EQUAL 0)
            message(WARNING "Failed to download Doxygen binary. Documentation target will not be available.")
            return()
        endif()

        # Extract the archive
        set(DOXYGEN_EXTRACT_DIR "${CMAKE_BINARY_DIR}/doxygen_bin")
        file(MAKE_DIRECTORY ${DOXYGEN_EXTRACT_DIR})

        if(APPLE)
            # Special handling for .dmg on macOS
            execute_process(
                COMMAND hdiutil attach ${DOXYGEN_ARCHIVE_FILE} -quiet -mountpoint ${DOXYGEN_EXTRACT_DIR}/mnt
                RESULT_VARIABLE RESULT
            )
            if(NOT RESULT EQUAL 0)
                message(WARNING "Failed to mount Doxygen .dmg. Documentation target will not be available.")
                return()
            endif()
            execute_process(
                COMMAND cp -R ${DOXYGEN_EXTRACT_DIR}/mnt/Doxygen.app ${DOXYGEN_EXTRACT_DIR}
                RESULT_VARIABLE RESULT
            )
            if(NOT RESULT EQUAL 0)
                message(WARNING "Failed to copy Doxygen.app. Documentation target will not be available.")
                return()
            endif()
            execute_process(
                COMMAND hdiutil detach ${DOXYGEN_EXTRACT_DIR}/mnt -quiet
            )
            set(DOXYGEN_EXECUTABLE "${DOXYGEN_EXTRACT_DIR}/Doxygen.app/Contents/MacOS/${DOXYGEN_BINARY_NAME}")
        else()
            # For zip/tar.gz on Windows/Linux
            file(ARCHIVE_EXTRACT INPUT ${DOXYGEN_ARCHIVE_FILE} DESTINATION ${DOXYGEN_EXTRACT_DIR})
            set(DOXYGEN_EXECUTABLE "${DOXYGEN_EXTRACT_DIR}/doxygen-${DOXYGEN_VERSION}/bin/${DOXYGEN_BINARY_NAME}")
        endif()

        if(NOT EXISTS ${DOXYGEN_EXECUTABLE})
            message(WARNING "Failed to find extracted Doxygen binary at ${DOXYGEN_EXECUTABLE}. Documentation target will not be available.")
            return()
        endif()

        message(STATUS "[setup_doxygen] Using pre-built Doxygen at ${DOXYGEN_EXECUTABLE}")
        set(DOXYGEN_FOUND TRUE)  # Mark as found
    else()
        message(STATUS "[setup_doxygen] Using system/installed Doxygen")
    endif()
endmacro()
