#ifndef CONFIG_H
#define CONFIG_H

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>

#include <vector>
#include <string>
#include <list>
#include <algorithm>

#include "../../options/user_options.h"
#include "../../utils/utils.h"

#include "../cost_function/cost_function.h"

// Parameters of the model
class Config
{
public:
    double l_d;                                     // Characteristic length {um}

    uint32_t num_terminals;                         // Current number of terminals
    uint32_t seed;                                  // Random seed for selecting points in the cloud       
    uint32_t np;                                    // Maximum number of feasible segments to be considered in the passive cost function evaluation
    uint32_t na;                                    // Maximum number of feasible segments to be considered in the active cost function evaluation

    double start_diameter;                          // Starting diameter for the segments
    double his_offset;                              // LAT offset from the His-Bundle

    double root_pos[3];                             // Root coordinates
    
    bool using_cloud_points;                        // Flag for using a cloud of points
    bool using_pmj_location;                        // Flag for using a PMJ location
    bool using_initial_network;                     // Flag for using an initial network

    std::string output_dir;                         // Location where the output files will be stored
    std::string cost_function_name;                 // Name of the cost function to be used
    std::string initial_network_name;               // Name of the iniitial network file
    std::string surface_filename;                   // Name of the file with surface triangles

    CostFunction *cost_fn;                          // Pointer to the cost function
    
    FILE *log_file;                                 // Reference to the LOGFILE

public:
    Config ();
    Config (User_Options *options);
    ~Config ();
    void print ();
private:
    void set_save_network (std::string output_dir);
    void set_parameters (User_Options *options);
    void set_save_network ();
    void set_cost_function (CostFunctionConfig *cost_fn_config);
};

#endif