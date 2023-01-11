set(zqPhysics_DIR ${CMAKE_CURRENT_LIST_DIR})
include(${CMAKE_CURRENT_LIST_DIR}/../Research-M/ResearchM/ResearchM.cmake)
get_filename_component(zqPhysics_INCLUDE_DIR ${CMAKE_CURRENT_LIST_DIR}/../ ABSOLUTE)
set(zqPhysics_INCLUDE_DIRS
	${zqPhysics_INCLUDE_DIR}
        ${ResearchM_INCLUDE_DIRS}
	${PROJECT_BINARY_DIR})

macro(zqPhysics_Deps project_name)
#	set(CMAKE_CXX_STANDARD 17)
#	link_directories(${zqPhysics_DIR}/../deps/lib/win64)	
#	if(NOT EXISTS ${zqPhysics_DIR}/../deps/inc)
#		message(FATAL_ERROR "deps not found, update submodule to get zqLibraryDeps")
#	else()
#		message(STATUS "zqLibraryDeps found as submodule")
#		include_directories(${zqPhysics_DIR}/../deps/inc)
#		if(${MSVC}) 
#			set(BIN_PATH "${zqPhysics_DIR}/../deps/bin/win64;${zqPhysics_DIR}/../deps/bin/win32")
#		endif()
#	endif()
#	# begin to find opengl
#	if(${MSVC})
#		# we already included OpenGL, so don't need to do anything more
#	else()
#		find_package(OpenGL REQUIRED)
#		find_package(GLUT REQUIRED)
#		find_package(GLEW REQUIRED)
#		target_link_libraries(${project_name} ${OPENGL_LIBRARIES} ${GLUT_LIBRARY})
#		target_link_libraries(${project_name} ${GLEW_LIBRARIES})
#	endif()
endmacro()	

	
