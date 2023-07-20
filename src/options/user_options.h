//
// Created by bergolho on 12/02/19.
//

#ifndef USER_OPTIONS_H
#define USER_OPTIONS_H

#include <cstdio>
#include <cstdlib>
#include <cstdbool>
#include <cstdint>
#include <cstring>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>

#include "../ini_parser/ini.h"
#include "../ini_parser/ini_file_sections.h"

#include "cost_function_config.h"
#include "pmj_config.h"
#include "cloud_config.h"

#define PRINT_LINE "====================================================================================="
#define PRINT_DOTS "....................................................................................."

class User_Options
{
public:
    uint32_t seed;
    
    double root_pos[3];

    double start_diameter;
    double his_offset;
    double characteristic_length;

    uint32_t np;
    uint32_t na;

    bool use_cloud_points;
    CloudConfig *cloud_config;

    bool use_pmj_location;
    PMJConfig *pmj_config;

    bool use_initial_network;
    std::string initial_network_filename;

    std::string output_dir;

    CostFunctionConfig *cost_function_config;

public:
    User_Options (const char filename[]);
    ~User_Options ();
    void read_config_file (const char filename[]);
    void print ();
};

int parse_config_file(void *user, const char *section, const char *name, const char *value);

#endif