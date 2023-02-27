//
// Created by bergolho on 20/09/21.
//

#ifndef GRAPH_H
#define GRAPH_H

#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <cstring>
#include <vector>
#include <queue>
#include <algorithm>

#include "../network/segment.h"

class GraphNode;
class GraphEdge;

class GraphNode
{
public:
    u_int32_t id;
    u_int32_t num_edges;
    double pos[3];
    GraphEdge *list_edges;
    GraphNode *next;
public:
    GraphNode (const u_int32_t id, const double pos[]);
    ~GraphNode ();
    void print ();
};

class GraphEdge
{
public:
    u_int32_t id;
    double length;
    char type;
    GraphNode *dest;
    GraphEdge *next;
public:
    GraphEdge (const u_int32_t id, const double length, GraphNode *dest);
    ~GraphEdge ();
    void print ();
};

class Graph
{
public:
    GraphNode *list_nodes;
    GraphNode *last_node;
    u_int32_t total_nodes;
    u_int32_t total_edges;
    double dx;
public:
    Graph ();
    ~Graph ();
    void insert_node (const double pos[]);
    void insert_edge (const uint32_t src_id, const uint32_t dst_id);
    GraphNode* search_node (const uint32_t id);
    void depth_first_search (GraphNode *u, const uint32_t level, uint32_t *map_network_to_mesh, bool *dfs_visited);
    void grow_segment (GraphNode *u, GraphEdge *v, uint32_t *map_network_to_mesh);
    void build_and_sort_segments (std::vector< std::pair<uint32_t,uint32_t> > &sorted_segments);
    void build_shocker_structure_from_graph(std::vector< std::pair<uint32_t,uint32_t> > sorted_segments,\
                                            const double diameter, const double lat_offset,\
                                            std::vector<Node*> &p_nodes, std::vector<Segment*> &p_segments);
    void print ();
    void write (std::string output_dir);
private:
    void build_shocker_nodes (std::vector<Node*> &nodes);
    void build_shocker_segments (std::vector< std::pair<uint32_t,uint32_t> > sorted_segments, const double diameter,\
                                std::vector<Node*> nodes, std::vector<Segment*> &segments);
    void update_diameter_and_LAT(std::vector<Segment*> &segments, const double diameter, const double lat_offset);
    void update_flags (std::vector<Node*> &nodes, std::vector<Segment*> &segments);
    void free_list_nodes (GraphNode *head);
    void free_list_edges (GraphNode *node);

};

Graph* convert_shocker_to_graph (std::vector<Node*> node_list, std::vector<Segment*> segment_list);
Graph* build_mesh_purkinje (Graph *network, const double dx);
void calc_unitary_vector (double d_ori[], GraphNode *u, GraphNode *v);

#endif