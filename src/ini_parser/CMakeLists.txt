set(SOURCE_FILES_INI_PARSER
        ini_file_sections.h
        ini.cpp
        ini.h
        )

# Add a static library for the linker of the main program
ADD_LIBRARY(ini_parser STATIC ${SOURCE_FILES_INI_PARSER})
