set(zqPhysics_DIR ${CMAKE_CURRENT_LIST_DIR})
include(${CMAKE_CURRENT_LIST_DIR}/../Research-M/ResearchM/ResearchM.cmake)
get_filename_component(zqPhysics_INCLUDE_DIR ${CMAKE_CURRENT_LIST_DIR}/../ ABSOLUTE)
set(zqPhysics_INCLUDE_DIRS
	${zqPhysics_INCLUDE_DIR}
    ${ResearchM_INCLUDE_DIRS}
	${PROJECT_BINARY_DIR})

	
