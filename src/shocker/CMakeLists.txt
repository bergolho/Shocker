set(SOURCE_FILES_SHOCKER shocker.h shocker.cpp
                         config/config.h config/config.cpp
                         cost_function/cost_function.h cost_function/cost_function.cpp
                         logger/logger.h logger/logger.cpp
                         error/error.h error/error.cpp
                         cloud/cloud.h cloud/cloud.cpp
                         cloud/point.h cloud/point.cpp
                         cloud/pmj.h cloud/pmj.cpp
                         network/node.h network/node.cpp
                         network/segment.h network/segment.cpp
                         graph/graph.h graph/graph.cpp
                        constants.h)

ADD_LIBRARY(shocker STATIC ${SOURCE_FILES_SHOCKER})

TARGET_LINK_LIBRARIES(shocker options utils ${VTK_LIBRARIES})

if(NOT MSVC)
    set_property(TARGET shocker PROPERTY POSITION_INDEPENDENT_CODE ON)
endif()
