cmake_minimum_required (VERSION 3.2)

message("This CMakeList.txt is written for 64bit system only")

project (ResearchPExamples LANGUAGES CXX)
set(ResearchP_ROOT_DIR ${CMAKE_SOURCE_DIR})


if(${MSVC}) 
	message(STATUS "Build with MSVC")
elseif(${UNIX})
	message(STATUS "Build with UNIX or LINUX")
else()
	message(WARNING "Build with unknown compiler")
endif()


###############################################
# Settings
###############################################
set_property(GLOBAL PROPERTY USE_FOLDERS ON)
add_definitions(-DUNICODE -D_UNICODE)

if(NOT CMAKE_BUILD_TYPE) 
	set(CMAKE_BUILD_TYPE Release)
endif() 
set(CMAKE_DEBUG_POSTFIX "d")

if(${MSVC}) 
	SET(CMAKE_ARCHIVE_OUTPUT_DIRECTORY_RELEASE  ${PROJECT_BINARY_DIR}/bin/Release)
	SET(CMAKE_ARCHIVE_OUTPUT_DIRECTORY_DEBUG  ${PROJECT_BINARY_DIR}/bin/Debug)
	SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE  ${PROJECT_BINARY_DIR}/bin/Release)
	SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG  ${PROJECT_BINARY_DIR}/bin/Debug)
endif()

###############################################
# Config data path
###############################################
set(MATH_DATA_PATH "${CMAKE_SOURCE_DIR}/Research-M/data")
configure_file(
	"${CMAKE_SOURCE_DIR}/Research-M/examples/math_data_path.h.in" 
	"${PROJECT_BINARY_DIR}/math_data_path.h"
	@ONLY)
message(STATUS "${PROJECT_BINARY_DIR}/data_path.h created,\n DATA_PATH=${MATH_DATA_PATH}")

###############################################
# Options and config
###############################################
include(ResearchP/ResearchP.cmake)
include_directories(${ResearchP_INCLUDE_DIRS})
message(STATUS "==============================")
# ---------- CUDA ----------
if(yzLib_ENABLE_CUDA)
	find_package(CUDA  REQUIRED)

	if(CUDA_FOUND)
		message(STATUS "CUDA version: ${CUDA_VERSION}")
		include_directories(${CUDA_INCLUDE_DIRS})  
		link_directories(${CUDA_TOOLKIT_ROOT_DIR}/lib/x64)
	else()
		message(ERROR "CUDA not found, you need to install CUDA, or disable yzLib_ENABLE_CUDA")
	endif()
endif()
###############################################
# Example Projects
###############################################
macro(SUBDIRLIST result curdir)
	file(GLOB children RELATIVE ${curdir} ${curdir}/*)
	set(dirlist "")
	foreach(child ${children})
	if(IS_DIRECTORY ${curdir}/${child})
		list(APPEND dirlist ${child})
	endif()
	endforeach()
	set(${result} ${dirlist})
endmacro()

SUBDIRLIST(example_groups ${CMAKE_SOURCE_DIR}/examples)

foreach(group ${example_groups})
	if(EXISTS ${CMAKE_SOURCE_DIR}/examples/${group}/CMakeLists.txt)
		set(example ${group})
		add_subdirectory(examples/${group})				
	else()
		SUBDIRLIST(examples ${CMAKE_SOURCE_DIR}/examples/${group})
		foreach(example ${examples})
			if(EXISTS ${CMAKE_SOURCE_DIR}/examples/${group}/${example}/CMakeLists.txt)
				add_subdirectory(examples/${group}/${example})		
			endif()
		endforeach()
	endif()
endforeach()
