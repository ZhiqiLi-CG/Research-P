set(ResearchP_DIR ${CMAKE_CURRENT_LIST_DIR})
include(${ResearchP_DIR}/../zqPhysics/zqPhysics.cmake)
include_directories(${zqPhysics_INCLUDE_DIRS})
get_filename_component(ResearchP_INCLUDE_DIR ${CMAKE_CURRENT_LIST_DIR}/../ ABSOLUTE)
set(ResearchP_INCLUDE_DIRS
	${ResearchP_INCLUDE_DIR}
	${PROJECT_BINARY_DIR})
	
