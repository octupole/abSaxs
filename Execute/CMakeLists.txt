# Policy
cmake_minimum_required( VERSION 3.0) 


if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()
set(SOURCE TrjRead.cpp  trjInput.cpp ExecabSaxs.cpp)
		
add_library(Execute STATIC ${SOURCE})
target_include_directories(Execute PUBLIC "${PROJECT_SOURCE_DIR}/abInitio" "${PROJECT_SOURCE_DIR}/utils" 
	"${PROJECT_SOURCE_DIR}/include"  "${PROJECT_SOURCE_DIR}/fftw3pp") 
	