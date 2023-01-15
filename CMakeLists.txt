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
set(PHYSICS_DATA_PATH "${CMAKE_SOURCE_DIR}/data")
configure_file(
	"${CMAKE_SOURCE_DIR}/examples/physics_data_path.h.in" 
	"${PROJECT_BINARY_DIR}/physics_data_path.h"
	@ONLY)
message(STATUS "${PROJECT_BINARY_DIR}/physics_path.h created,\n DATA_PATH=${PHYSICS_DATA_PATH}")

###############################################
# Options and config
###############################################
include(ResearchP/ResearchP.cmake)
include_directories(${ResearchP_INCLUDE_DIRS})
Set_ResearchP_Env()

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
