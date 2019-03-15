################################################################################
include(DownloadProject)

# With CMake 3.8 and above, we can hide warnings about git being in a
# detached head by passing an extra GIT_CONFIG option
if(NOT (${CMAKE_VERSION} VERSION_LESS "3.8.0"))
    set(TMESH_KERNEL_EXTRA_OPTIONS "GIT_CONFIG advice.detachedHead=false")
else()
    set(TMESH_KERNEL_EXTRA_OPTIONS "")
endif()

# Shortcut function
function(tmesh_kernel_download_project name)
    download_project(
        PROJ         ${name}
        SOURCE_DIR   ${TMESH_KERNEL_EXTERNAL}/${name}
        DOWNLOAD_DIR ${TMESH_KERNEL_EXTERNAL}/.cache/${name}
        QUIET
        ${TMESH_KERNEL_EXTRA_OPTIONS}
        ${ARGN}
    )
endfunction()

################################################################################

## Catch2 BSL 1.0
function(tmesh_kernel_download_catch2)
    tmesh_kernel_download_project(Catch2
        GIT_REPOSITORY https://github.com/catchorg/Catch2.git
        GIT_TAG        v2.4.2
    )
endfunction()

## tbb Apache-2.0
function(tmesh_kernel_download_tbb)
    tmesh_kernel_download_project(tbb
        GIT_REPOSITORY https://github.com/wjakob/tbb.git
        GIT_TAG        08b4341a1893a72656467e96137f1f99d0112547
    )
endfunction()


## Sanitizers MIT
function(tmesh_kernel_download_sanitizers)
    tmesh_kernel_download_project(sanitizers-cmake
        GIT_REPOSITORY https://github.com/arsenm/sanitizers-cmake.git
        GIT_TAG        6947cff3a9c9305eb9c16135dd81da3feb4bf87f
    )
endfunction()
