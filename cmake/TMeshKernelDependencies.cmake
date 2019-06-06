# Prepare dependencies
#
# For each third-party library, if the appropriate target doesn't exist yet,
# download it via external project, and add_subdirectory to build it alongside
# this project.

### Configuration
set(TMESH_KERNEL_ROOT     "${CMAKE_CURRENT_LIST_DIR}/..")
set(TMESH_KERNEL_EXTERNAL "${TMESH_KERNEL_ROOT}/external")

# Download and update 3rdparty libraries
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR})
list(REMOVE_DUPLICATES CMAKE_MODULE_PATH)
include(TMeshKernelDownloadExternal)

################################################################################
# Required libraries
################################################################################

# Sanitizers
if(TMESH_KERNEL_WITH_SANITIZERS)
    tmesh_kernel_download_sanitizers()
    find_package(Sanitizers)
endif()

tmesh_kernel_download_catch2()
