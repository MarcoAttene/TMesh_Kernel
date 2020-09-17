################################################################################
cmake_minimum_required(VERSION 3.1)
project(TMesh_Kernel)
################################################################################

# Uncomment and adapt the following definition after having downloaded
# mpir and compiled it for 64bit.
#set(MPIR_LOCATION "D:/SOFTWARE/mpir-master/msvc/vs19/lib_mpir_gc/x64/Release")

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
option(USE_LAZY_KERNEL   "Use lazy kernel"       ON)

################################################################################
# library
################################################################################

if(${USE_LAZY_KERNEL})
	if(NOT DEFINED MPIR_LOCATION)
	  message(FATAL_ERROR "MPIR_LOCATION is undefined. Please set it on line 8 in CMakeLists.txt")
	endif(NOT DEFINED MPIR_LOCATION)
	set(TMESH_LIB_PATH "kernel_Lazy64")
else()
	set(TMESH_LIB_PATH "kernel_Fast64")
endif()

add_library(${TMESH_LIB_PATH} STATIC
	src/coordinates.cpp
	src/expansion.cpp
	src/mixedPredicates.cpp
	src/orientation.cpp
	src/point.cpp
	src/tmesh.cpp
	)

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
	target_compile_options(${TMESH_LIB_PATH} PUBLIC -frounding-math)
	set(MPIR_LIB ${MPIR_LOCATION}/mpir.a)
elseif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	target_compile_options(${TMESH_LIB_PATH} PUBLIC "/fp:strict")
	target_compile_definitions(${TMESH_LIB_PATH} PUBLIC -D_CRT_SECURE_NO_WARNINGS)
	set(MPIR_LIB ${MPIR_LOCATION}/mpir.lib)
endif()

if(${USE_LAZY_KERNEL})
	target_compile_definitions(${TMESH_LIB_PATH} PUBLIC -DUSE_LAZY_KERNEL)
	target_compile_definitions(${TMESH_LIB_PATH} PUBLIC -DUSE_HYBRID_KERNEL)
endif()

if(CMAKE_SIZEOF_VOID_P EQUAL 8)
	target_compile_definitions(${TMESH_LIB_PATH} PUBLIC -DIS64BITPLATFORM)
endif()

# Public include directory
target_include_directories(${TMESH_LIB_PATH} PUBLIC include)

# Use C++14
target_compile_features(${TMESH_LIB_PATH} PUBLIC ${CXX14_FEATURES})


################################################################################
# Required libraries
################################################################################

if(${USE_LAZY_KERNEL})
	target_include_directories(${TMESH_LIB_PATH} PUBLIC ${MPIR_LOCATION}/)
	target_link_libraries(${TMESH_LIB_PATH} PUBLIC ${MPIR_LIB})
endif()
