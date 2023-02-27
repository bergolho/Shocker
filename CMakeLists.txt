cmake_minimum_required(VERSION 3.0.1)
project(Shocker)

SET(CMAKE_BUILD_TYPE "Debug")
SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)
SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/shared-libs)

FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${VTK_USE_FILE})

ADD_SUBDIRECTORY(src/options)
ADD_SUBDIRECTORY(src/ini_parser)
ADD_SUBDIRECTORY(src/utils)
ADD_SUBDIRECTORY(src/shocker)

ADD_EXECUTABLE(Shocker src/main.cpp)

ADD_DEPENDENCIES(Shocker shocker options ini_parser utils)

TARGET_LINK_LIBRARIES(Shocker ${VTK_LIBRARIES} shocker options ini_parser utils dl m)
