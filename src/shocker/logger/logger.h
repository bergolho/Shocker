#ifndef LOGGER_H
#define LOGGER_H

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

class Logger
{
public:
    uint32_t iter;
public:
    Logger (User_Options *options);
    ~Logger ();
    void write_lmin_to_file (const std::string output_dir, const uint32_t num_terminal, const double l_min);
    void write_full_network_to_vtk (const std::vector<Node*> &n_list, const std::vector<Segment*> &s_list,\
                                    const std::string output_dir, const double his_offset, const uint32_t num_terminals);
    void write_full_network_monoalg_to_vtk (const std::vector<Node*> &n_list, const std::vector<Segment*> &s_list,\
                                    const std::string output_dir, const double his_offset, const uint32_t num_terminals);
    void write_minimum_network_to_vtk (const std::vector<Node*> &n_list, const std::vector<Segment*> &s_list,\
                                    const std::string output_dir, const double his_offset, const uint32_t num_terminals);
    void write_minimum_network_monoalg_to_vtk (const std::vector<Node*> &n_list, const std::vector<Segment*> &s_list,\
                                    const std::string output_dir, const double his_offset, const uint32_t num_terminals);
    void write_simulation_time (const std::string output_dir, const long res_time, const uint32_t inactive_points_time, const uint32_t active_points_time,\
                                const uint32_t inactive_search_time, const uint32_t active_search_time,\
                                const uint32_t inactive_eval_time, const uint32_t active_eval_time,\
                                const uint32_t main_loop_time, const uint32_t force_connection_no_lat_error_tolerance_time, const uint32_t force_connection_no_distance_criterion_time,\
                                const uint32_t force_connection_no_distance_criterion_geodesic_time, const uint32_t force_connection_no_distance_criterion_line_time,\
                                const uint32_t geodesic_pathway_time,\
                                const double conv_rate);
    void write_configuration_file (User_Options *options);
};

#endif