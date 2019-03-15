################################################################################
# Find MPIR and build it as part of the current build
################################################################################

if(TARGET mpir)
	return()
endif()

################################################################################


add_library(mpir INTERFACE)

if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	target_include_directories(mpir INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/mpir)

	if(${TMESH_KERNEL_IS_Pentium_4_SSE2})
		target_link_libraries(mpir INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/mpir/mpir32.lib)
	elseif(${TMESH_KERNEL_IS_intel_core2})
		target_link_libraries(mpir INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/mpir/mpir64.lib)
	elseif(CMAKE_SIZEOF_VOID_P EQUAL 8)
		target_link_libraries(mpir INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/mpir/mpir64-generic.lib)
	else()
		target_link_libraries(mpir INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/mpir/mpir32-generic.lib)
	endif()
elseif(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
	if(${TMESH_KERNEL_USE_MPIR_PROVIDED_BINARY})
		target_include_directories(mpir INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/mpir)
		target_link_libraries(mpir INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/mpir/libmpir.a)
	else()
		find_package(MPIR)

		if(NOT ${MPIR_FOUND})
			MESSAGE(FATAL_ERROR "Unable to find mpir")
		endif()

		target_include_directories(mpir INTERFACE ${MPIR_INCLUDE_DIR})
		target_link_libraries(mpir INTERFACE ${MPIR_LIBRARIES})
	endif()
elseif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	find_package(MPIR)

	if(NOT ${MPIR_FOUND})
		MESSAGE(FATAL_ERROR "Unable to find mpir")
	endif()

	target_include_directories(mpir INTERFACE ${MPIR_INCLUDE_DIR})
	target_link_libraries(mpir INTERFACE ${MPIR_LIBRARIES})
endif()

