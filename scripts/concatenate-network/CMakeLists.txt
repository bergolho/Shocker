cmake_minimum_required(VERSION 3.0.1)
project(ConcatenateNetwork)

if (CMAKE_VERSION VERSION_LESS "3.1")
    if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
      set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++11")
    endif ()
else ()
    set (CMAKE_CXX_STANDARD 11)
endif ()

SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)

SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -fno-strict-aliasing -Wall -Wno-unused-function")

ADD_EXECUTABLE(ConcatenateNetwork src/main.cpp)

TARGET_LINK_LIBRARIES(ConcatenateNetwork ${VTK_LIBRARIES} m)