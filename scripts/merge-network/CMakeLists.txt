cmake_minimum_required(VERSION 3.0.1)
project(MergeNetwork)

SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)
SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/shared-libs)

ADD_SUBDIRECTORY(src/reader)
ADD_SUBDIRECTORY(src/graph)
ADD_SUBDIRECTORY(src/utils)

ADD_EXECUTABLE(MergeNetwork src/main.cpp)

ADD_DEPENDENCIES(MergeNetwork reader graph utils)

TARGET_LINK_LIBRARIES(MergeNetwork reader graph utils dl m)