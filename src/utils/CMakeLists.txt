set(SOURCE_FILES_UTILS utils.h utils.cpp
    stop_watch.cpp stop_watch.h)

ADD_LIBRARY(utils STATIC ${SOURCE_FILES_UTILS})

if(NOT MSVC)
    set_property(TARGET utils PROPERTY POSITION_INDEPENDENT_CODE ON)
endif()