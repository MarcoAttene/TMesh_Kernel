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
    polyfem_download_sanitizers()
    find_package(Sanitizers)
endif()

# TBB
if(${TMESH_KERNEL_ENABLE_TBB} AND NOT TARGET tbb_static)
	tmesh_kernel_download_tbb()
	set(TBB_BUILD_STATIC ON CACHE BOOL " " FORCE)
	set(TBB_BUILD_SHARED OFF CACHE BOOL " " FORCE)
	set(TBB_BUILD_TBBMALLOC OFF CACHE BOOL " " FORCE)
	set(TBB_BUILD_TBBMALLOC_PROXY OFF CACHE BOOL " " FORCE)
	set(TBB_BUILD_TESTS OFF CACHE BOOL " " FORCE)

	add_subdirectory(${THIRD_PARTY_DIR}/tbb tbb)
	set_property(TARGET tbb_static tbb_def_files PROPERTY FOLDER "dependencies")
	target_compile_definitions(tbb_static PUBLIC -DTMESH_KERNEL_USE_TBB)
endif()

tmesh_kernel_download_catch2()
