#ifndef ERROR_H
#define ERROR_H

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>

#include <iostream>
#include <vector>
#include <string>
#include <list>
#include <algorithm>

#include "../../options/user_options.h"
#include "../../utils/utils.h"

#include "../network/segment.h"
#include "../cloud/cloud.h"

class Error
{
public:
    uint32_t counter_main_loop_connections;                         // Number of PMJs connected in the main loop (Default)
    uint32_t counter_no_lat_error_tolerance_connections;            // Number of PMJs connected after droping the LAT error tolerance (Pos-Processing)
    uint32_t counter_no_distance_criterion_connections;             // Number of PMJs connected after droping the distance criterion  (Pos-Processing)
    uint32_t counter_no_distance_criterion_connections_geodesic;    // Number of PMJs connected after droping the distance criterion and connected with geodesic pathway  (Pos-Processing)
    uint32_t counter_no_distance_criterion_connections_line;        // Number of PMJs connected after droping the distance criterion and connected with straight line  (Pos-Processing)
    double epsilon_2ms;                                     // Percentage of PMJs that are connected with an error below 2ms
    double epsilon_5ms;                                     // Percentage of PMJs that are connected with an error below 5ms
    double rmse;                                            // Root Mean Square Error of the connected PMJs
    double rrmse;                                           // Relative Root Mean Square Error of the connected PMJs
    double max_lat_error;                                   // Maximum LAT error of the connected PMJs
    double min_max_aprox_lat[2];                            // Minimum/Maximum LAT of the generated network
    double min_max_ref_lat[2];                              // Minimum/Maximum LAT of the reference network
    double min_max_term_lat[2];                             // Minimum/Maximum LAT of the connected PMJs in the generated network

public:
    Error ();
    ~Error ();
    Error* copy ();
    void update_min_max_terminal_lat (std::vector<Segment*> s_list);
    void calculate_electric_error (PMJ_Data *pmj_data);
    void concatenate (Error *input);
    void reset ();
    void print ();
private:
};

#endif