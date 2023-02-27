set(SOURCE_FILES_SEGMENT segment.h segment.cpp
                         point.h point.cpp)

ADD_LIBRARY(network STATIC ${SOURCE_FILES_SEGMENT})

TARGET_LINK_LIBRARIES(network utils)

if(NOT MSVC)
    set_property(TARGET network PROPERTY POSITION_INDEPENDENT_CODE ON)
endif()