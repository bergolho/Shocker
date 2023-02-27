cmake_minimum_required(VERSION 3.0.1)
project(NetworkReader)

SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)
SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/shared-libs)

ADD_SUBDIRECTORY(src/reader)
ADD_SUBDIRECTORY(src/graph)
ADD_SUBDIRECTORY(src/network)
ADD_SUBDIRECTORY(src/utils)

ADD_EXECUTABLE(NetworkReader src/main.cpp)

ADD_DEPENDENCIES(NetworkReader reader graph network utils)

TARGET_LINK_LIBRARIES(NetworkReader reader graph network utils dl m)
