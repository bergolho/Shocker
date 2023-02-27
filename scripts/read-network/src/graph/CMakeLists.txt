set(SOURCE_FILES_OPTIONS graph.h graph.cpp
                         pqueue.h pqueue.cpp
                         constants.h)

ADD_LIBRARY(graph STATIC ${SOURCE_FILES_OPTIONS})

if(NOT MSVC)
    set_property(TARGET graph PROPERTY POSITION_INDEPENDENT_CODE ON)
endif()