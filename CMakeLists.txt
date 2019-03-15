################################################################################
cmake_minimum_required(VERSION 3.1)
project(TMesh_Kernel)
################################################################################

# cmake_policy(SET CMP0063 NEW)
# set(CMAKE_CXX_VISIBILITY_PRESET hidden)

# Detects whether this is a top-level project
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
	set(TMESH_KERNEL_TOPLEVEL_PROJECT ON)
else()
	set(TMESH_KERNEL_TOPLEVEL_PROJECT OFF)
endif()

set(THIRD_PARTY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
list(APPEND CMAKE_MODULE_PATH ${THIRD_PARTY_DIR}/Catch2/contrib)

# Color output
include(UseColors)


# Extra warnings
include(Warnings)

# Use C++11/14
include(CXXFeatures)

# Sort projects inside the solution
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

# Generate position independent code by default
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

################################################################################

if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	foreach(config ${CMAKE_CONFIGURATION_TYPES})
		string(TOUPPER ${config} config)
		string(REPLACE /MD /MT CMAKE_C_FLAGS_${config} "${CMAKE_C_FLAGS_${config}}")
		string(REPLACE /MD /MT CMAKE_CXX_FLAGS_${config} "${CMAKE_CXX_FLAGS_${config}}")
	endforeach()
endif()

################################################################################


# Options for enabling/disabling optional libraries
option(TMESH_KERNEL_ENABLE_TBB         "Enable TBB"             ON)

# Sanitizer options
option(TMESH_KERNEL_SANITIZE_ADDRESS   "Sanitize Address"       OFF)
option(TMESH_KERNEL_SANITIZE_MEMORY    "Sanitize Memory"        OFF)
option(TMESH_KERNEL_SANITIZE_THREAD    "Sanitize Thread"        OFF)
option(TMESH_KERNEL_SANITIZE_UNDEFINED "Sanitize Undefined"     OFF)


option(TMESH_KERNEL_USE_HYBRID_KERNEL "Use hybrid kernel"     ON)
option(TMESH_KERNEL_USE_LAZY_KERNEL   "Use lazy kernel"       ON)
option(TMESH_KERNEL_EXTENSIBLE_TMESH  "Enable polymorphism"   OFF)

if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	option(TMESH_KERNEL_IS_Pentium_4_SSE2  "Uses MPIR optimized for processors of the Pentium 4 SSE2 family"         OFF)
	option(TMESH_KERNEL_IS_intel_core2     "Uses MPIR optimized for processors of the Intel(TM) core-2(TM) family"   OFF)
elseif(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
	option(TMESH_KERNEL_USE_MPIR_PROVIDED_BINARY     "Uses the MPIR precompiled on a generic 64 bit Ubuntu Linux"   ON)
endif()


################################################################################
# Dependencies
################################################################################

# Sanitizers
if(TMESH_KERNEL_WITH_SANITIZERS)
	list(APPEND CMAKE_MODULE_PATH ${THIRD_PARTY_DIR}/sanitizers-cmake/cmake)
endif()

# Setup dependencies
include(TMeshKernelDependencies)

################################################################################
# library
################################################################################

add_library(${PROJECT_NAME}
	src/coordinates.cpp
	src/expansion.cpp
	src/mixedPredicates.cpp
	src/orientation.cpp
	src/point.cpp
	src/tmesh.cpp
	)

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
	target_compile_options(${PROJECT_NAME} PUBLIC -frounding-math)
elseif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	target_compile_options(${PROJECT_NAME} PUBLIC "/fp:strict")
endif()

if(${TMESH_KERNEL_USE_HYBRID_KERNEL})
	target_compile_definitions(${PROJECT_NAME} PUBLIC -DUSE_HYBRID_KERNEL)
endif()

if(${TMESH_KERNEL_USE_LAZY_KERNEL})
	target_compile_definitions(${PROJECT_NAME} PUBLIC -DUSE_LAZY_KERNEL)
	target_compile_definitions(${PROJECT_NAME} PUBLIC -DUSE_HYBRID_KERNEL)
endif()

if(CMAKE_SIZEOF_VOID_P EQUAL 8)
	target_compile_definitions(${PROJECT_NAME} PUBLIC -DIS64BITPLATFORM)
endif()

# Public include directory
target_include_directories(${PROJECT_NAME} PUBLIC include)

# Extra warnings
target_link_libraries(${PROJECT_NAME} PRIVATE warnings::all)

# Use C++14
target_compile_features(${PROJECT_NAME} PUBLIC ${CXX14_FEATURES})


if(TMESH_KERNEL_WITH_SANITIZERS)
	add_sanitizers(${PROJECT_NAME})
endif()


################################################################################
# Required libraries
################################################################################

# TBB
if(TMESH_KERNEL_ENABLE_TBB)
	target_include_directories(${PROJECT_NAME} SYSTEM PUBLIC ${THIRD_PARTY_DIR}/tbb/include)
	target_link_libraries(${PROJECT_NAME} PUBLIC tbb_static)
endif()


# Catch2
add_library(catch INTERFACE)
target_include_directories(catch SYSTEM INTERFACE ${THIRD_PARTY_DIR}/Catch2/single_include/catch2)


include(mpir)
target_link_libraries(${PROJECT_NAME} PUBLIC mpir)

################################################################################
# Subdirectories
################################################################################

# Compile extras only if this is a top-level project
if(TMESH_KERNEL_TOPLEVEL_PROJECT)
	# Unit tests
	include(CTest)
	enable_testing()
	add_subdirectory(tests)
endif()
