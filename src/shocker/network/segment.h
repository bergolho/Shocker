//
// Created by bergolho on 18/08/21.
//

#ifndef SEGMENT_H
#define SEGMENT_H

#include <cstdio>
#include <cstdlib>
#include <cstdbool>
#include <cstdint>
#include <cmath>

#include "../constants.h"
#include "../../utils/utils.h"
#include "../../utils/stop_watch.h"
#include "node.h"

class Segment
{
public:
    uint32_t id;
    double diameter;
    double length;
    double lat;
    double middle_pos[3];
    
    bool mintree;
    
    Node *src;
    Node *dest;

    Segment *left;
    Segment *right;
    Segment *parent;
public:
    Segment ();
    Segment (const uint32_t id, const double diameter, const double lat, const bool mintree,\
            Node *src, Node *dest,\
            Segment *left, Segment *right, Segment *parent);
    Segment (Segment *in);
    ~Segment ();
    void calc_middle_point (double pos[]);
    void calc_unitary_vector (double u[]);
    void tag_pathway ();
    double calc_radius ();
    double calc_length ();
    double calc_propagation_velocity ();
    double calc_local_activation_time ();
    double calc_terminal_local_activation_time ();
    double calc_pathway_length ();
    double calc_pathway_length (std::vector<double> &out);
    double calc_diameter (const double cv);
    double calc_conductivity ();
    double calc_dproj (const double pos[]);
    double calc_dortho (const double pos[]);
    double calc_dend (const double pos[]);
    uint32_t calc_level ();
    bool is_segment ();
    bool is_terminal ();
    bool is_bifurcation ();
    bool is_inside_region (const double center[], const double radius);
    bool is_mintree ();
    void print ();
};

Node* generate_bifurcation_node (const uint32_t num_nodes, Segment *iconn);
Node* generate_terminal_node (const uint32_t num_nodes, const double pos[]);
void move_bifurcation_position (Segment *iconn, Segment *ibiff, Segment *inew, const double pos[]);
void update_segment_pointers (Segment *iconn, Segment *ibiff, Segment *inew, const bool is_restore);
void recalculate_length (std::vector<Segment*> &s_list);
void recalculate_local_activation_time (std::vector<Segment*> &s_list);
void recalculate_length_branch (Segment *inew, Segment *ibiff);
void recalculate_local_activation_time_branch (Segment *inew, Segment *ibiff);
void get_segment_length (std::vector<Segment*> s_list, std::vector<double> &segments);
void get_bifurcation_angles(std::vector<Segment*> s_list, std::vector<double> &angles);
void compute_distance(std::vector<Segment*> s_list, const double pos[], std::vector< std::pair<double,uint32_t> > &dist_array);
Segment* get_closest_segment_by_distance (std::vector<Segment*> s_list, const double pos[]);
std::vector<Segment*> get_closest_segment_list_by_distance (std::vector<Segment*> s_list, const double pos[], const uint32_t num_size);
std::vector<Segment*> get_closest_segment_list_by_LAT_error (std::vector<Segment*> s_list, const double pos[], const double ref_lat, const uint32_t num_size);
Segment* get_closest_segment_by_LAT (std::vector<Segment*> s_list, const double ref_lat);

void eliminate_segment_from_list (std::vector<Segment*> &s_list, Segment *s);
void order_segment_list (std::vector<Segment*> &s_list);
void print_segment_list (std::vector<Segment*> s_list);
void write_segment_list (const char filename[], std::vector<Segment*> s_list);
//void copy_segment_list (std::vector<Node*> in_points, std::vector<Segment*> in, std::vector<Segment*> &out);
//void concatenate_segment_list (std::vector<Segment*> in, std::vector<Node*> out_points, std::vector<Segment*> &out,\
                                const uint32_t offset_points, const uint32_t network_id);
                                
uint32_t get_total_time_local_activation_terminal ();

#endif