set(SOURCE_FILES_NETWORK segment.h segment.cpp
                         node.h node.cpp)

ADD_LIBRARY(network STATIC ${SOURCE_FILES_NETWORK})

TARGET_LINK_LIBRARIES(network utils)

if(NOT MSVC)
    set_property(TARGET network PROPERTY POSITION_INDEPENDENT_CODE ON)
endif()