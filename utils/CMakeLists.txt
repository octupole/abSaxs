# Policy
cmake_minimum_required( VERSION 3.0) 

if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()

set(SOURCE SGSmooth.cpp)

add_library(utils STATIC ${SOURCE})
target_include_directories(utils PUBLIC "${PROJECT_SOURCE_DIR}/include" "${PROJECT_SOURCE_DIR}/alglib") 

