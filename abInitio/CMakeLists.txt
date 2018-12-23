# Policy
cmake_minimum_required( VERSION 3.0) 

# Project name
if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()

set(SOURCE ABSaxs.cpp abRhoSaxs.cpp SaxsData.cpp Funktionell.cpp)
add_library(absaxs STATIC ${SOURCE})

target_include_directories(absaxs PUBLIC "${PROJECT_SOURCE_DIR}/Execute" "${PROJECT_SOURCE_DIR}/utils" 
	"${PROJECT_SOURCE_DIR}/include" "${PROJECT_SOURCE_DIR}/libalg" "${PROJECT_SOURCE_DIR}/fftw3pp"
	"${PROJECT_SOURCE_DIR}/libtraj" "${PROJECT_SOURCE_DIR}/Saxs") 

