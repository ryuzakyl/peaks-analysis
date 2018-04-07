cmake_minimum_required(VERSION 3.5)
project(peaks-analysis)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES src/config.hpp src/peaks-analysis.h src/peaks-analysis.cpp)
add_library(_peaks_analysis SHARED ${SOURCE_FILES})
