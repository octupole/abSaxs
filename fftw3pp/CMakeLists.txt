# Policy
cmake_minimum_required( VERSION 3.0) 


if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()

set(SOURCE Pfftw.cpp Prcfft3d.cpp Pcrfft3d.cpp)

add_library(fftw3pp STATIC ${SOURCE})
