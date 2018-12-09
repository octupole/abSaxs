# Policy
cmake_minimum_required( VERSION 3.0) 

# Project name
if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()

set(SOURCE ABSaxs.cpp RhoSaxs.cpp SaxsData.cpp)

add_library(absaxs STATIC ${SOURCE})
target_include_directories(absaxs PUBLIC "${PROJECT_SOURCE_DIR}/Execute" "${PROJECT_SOURCE_DIR}/utils" 
	"${PROJECT_SOURCE_DIR}/include""${PROJECT_SOURCE_DIR}/libalg") 

