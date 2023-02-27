//
// Created by bergolho on 17/08/21.
//

#ifndef SHOCKER_H
#define SHOCKER_H

#include <cstdio>
#include <cstdlib>
#include <cstdbool>
#include <cstdint>
#include <cmath>

#include <iostream>
#include <vector>
#include <string>
#include <algorithm>

#include "../options/user_options.h"
#include "../utils/utils.h"
#include "../utils/stop_watch.h"

#include "network/segment.h"
#include "config/config.h"
#include "cloud/cloud.h"
#include "error/error.h"
#include "graph/graph.h"
#include "logger/logger.h"

class Shocker
{
private:
    Config *params;
    Cloud *cloud;
    Error *error;
    Logger *logger;

    std::vector<Node*> node_list;
    std::vector<Segment*> segment_list;

public:
    Shocker(User_Options *options);
    ~Shocker();
    void grow ();
    bool generate_terminal ();
    bool generate_pmj (PMJPoint *pmj);
    bool attempt_pmj_connection ();
    void save_network_state ();
    void load_network_state ();
    void reset_network_state ();
    void print_network_info ();
    void write_full_network_to_vtk();
    void write_minimum_network_to_vtk ();
    void write_network_info_to_logfile (const double mean_segment_length, const double std_segment_length, const uint32_t ns,\
                                             const double mean_biff_angle, const double std_biff_angle, const uint32_t nb);
    std::vector<std::vector<double>> get_nodes_coordinates ();
private:
    void root_placement ();
    void make_root_default ();
    void grow_tree_using_cloud_points ();
    void generate_root_segment (const double x_prox[], const double x_dist[]);
    void update_and_check_integrity ();
    void build_minimum_network_maps (std::map<uint32_t,uint32_t> &nodes_map, std::map<uint32_t,uint32_t> &segments_map);
    void build_minimum_network_topology (std::map<uint32_t,uint32_t> nodes_map, std::map<uint32_t,uint32_t> segments_map,\
                                std::vector<Node> &min_network_nodes, std::vector<Segment> &min_network_segments);
    bool connection_search (const double pos[], const double l_min);
    bool fill_feasible_segments (std::vector<Segment*> &feasible_segments, const double pos[]);
    bool fill_feasible_segments_pmj (std::vector<Segment*> &feasible_segments, const double pos[], const double ref_lat);
    bool evaluate_cost_function_cloud_point (std::vector<Segment*> &feasible_segments, CloudPoint *point, Evaluation &best);
    bool evaluate_cost_function_pmj (std::vector<Segment*> &feasible_segments, PMJPoint *pmj, Evaluation &best);
    bool evaluate_pmj_local_activation_time (Segment *inew, PMJPoint *pmj);
    bool distance_criterion (Segment *s, const double pos[], const double l_min);
    bool process_pmjs (uint32_t &prev_num_connected_pmjs, uint32_t &cur_num_connected_pmjs);
    bool prune_inactive_segments ();
    bool check_minimum_distance (const double l_min);
    bool check_duplicates (std::vector<Node*> n_list, std::vector<CloudPoint> geo_points);
    bool pos_process ();
    bool force_connection_without_distance_criterion ();
    bool force_connection_without_distance_criterion_2 ();
    bool force_connection_without_lat_error_tolerance ();
    bool fill_feasible_pmjs (std::vector<uint32_t> &feasible_pmjs);
    bool attempt_diameter_adjustment ();
    uint32_t tag_minimum_network_segments ();
    double update_minimum_distance (uint32_t &tosses, const double l_min);
    Segment* build_branch (Segment *iconn, const double pos[], double &branch_size);
    Segment* prune_branch (Segment *inew);
    Segment* search_active_terminal_segment (PMJPoint *pmj);
    void recalculate_errors ();
    void extract_connected_pmjs ();
    bool refinement_attempt_diameter_adjustment ();
};

#endif