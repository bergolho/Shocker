#ifndef GRAPH_H
#define GRAPH_H

#include <cstdio>
#include <cstdlib>
#include <cstdbool>
#include <cstdint>
#include <cstring>
#include <vector>
#include <queue>
#include <map>
#include <string>

#include "../utils/utils.h"
#include "../reader/reader.h"
#include "../network/segment.h"

#include "pqueue.h"
#include "constants.h"

// --------------------------------------------------------------------------------
// 'PQUEUE' USER DEFINED STRUCTURES AND FUNCTIONS
typedef struct node_t
{
	pqueue_pri_t pri;
	uint32_t    val;
	size_t pos;
} node_t;


static int
cmp_pri(pqueue_pri_t next, pqueue_pri_t curr)
{
	//return (next < curr);         // equivalent to std::less<int>()
    return (next > curr);           // equivalent to std::greater<int>()
}


static pqueue_pri_t
get_pri(void *a)
{
	return ((node_t *) a)->pri;
}


static void
set_pri(void *a, pqueue_pri_t pri)
{
	((node_t *) a)->pri = pri;
}


static size_t
get_pos(void *a)
{
	return ((node_t *) a)->pos;
}


static void
set_pos(void *a, size_t pos)
{
	((node_t *) a)->pos = pos;
}
// --------------------------------------------------------------------------------

class Edge
{
public:
    bool is_minimum_tree;
    uint32_t dest_id;
    uint32_t branch_id;
    double length;
    double _lat;                // Parameter from Shocker VTK network file
    double _diameter;           // Parameter from Shocker VTK network file
    double _length;             // Parameter from Shocker VTK network file
    double _min_tree;           // Parameter from Shocker VTK network file
public:
    Edge (const uint32_t dest_id, const double length)
    {
        this->is_minimum_tree = false;
        this->dest_id = dest_id;
        this->length = length;
        this->branch_id = 0;
        this->_lat = 0.0;
        this->_diameter = 68.8608;
        this->_min_tree = 1.0;
    }
};

class Node
{
public:
    uint32_t id;
    double x, y, z;
    double lat;
    std::vector<Edge> list_edges;
public:
    Node () { }
    Node (const uint32_t id, const double x,const double y, const double z)
    {
        this->id = id;
        this->x = x;
        this->y = y;
        this->z = z;
    }
    void print ()
    {
        printf("|| [%u] (%g %g %g) {%g} ||",this->id,this->x,this->y,this->z,this->lat);
        for (uint32_t i = 0; i < this->list_edges.size(); i++)
            printf(" --> || [%u] %g %d [Branch = %u] || ",this->list_edges[i].dest_id,this->list_edges[i].length,(int)this->list_edges[i].is_minimum_tree,this->list_edges[i].branch_id);
        printf("\n");
    }
    bool is_terminal () { return (this->list_edges.size() == 0); }
};

class Graph
{
public:
    uint32_t root_index;                            // Index where the root point is located
    uint32_t total_nodes;                           // Total number of nodes
    uint32_t total_edges;                           // Total number of edges
    std::vector<Node> list_nodes;                   // List of nodes
    std::vector<double> dist;                       // Shortest distance of each cell from the root point
    std::vector<double> lat;                        // Local activation time of each cell
    std::vector<int> parent;                        // Parent indexes of each cell
    std::vector<int> dfs_counter;                   // DFS visit order counter
    std::vector<bool> dfs_visited;                  // DFS visit flag
    std::vector<uint32_t> terminals_indexes;        // Active PMJ's indexes
    double min_LAT;
public:
    Graph ();
    Graph (std::vector<Point> points, std::vector<Line> lines, std::vector<std::vector<double>> cell_data, std::string root_filename);
    uint32_t count_active_edges ();
    uint32_t get_closest_terminal_point (Point p, const bool is_reference);
    void depth_first_search_default (const uint32_t src_id);
    void depth_first_search_branches (const uint32_t src_id, std::vector<double> &the_branches);
    void depth_first_search_topology (const uint32_t src_id);
    void dfs_default (Node u, std::vector<bool> &dfs_visited);
    void dfs_branches (Node u, std::vector<double> &segments);
    void dfs_topology (Node *current);
    void dijkstra (const uint32_t src_id);
    void dijkstra_2 (const uint32_t src_id);
    void fill_terminal_indexes (const uint32_t root_index);
    void tag_minimum_network (std::vector<Point> pmjs);
    uint32_t tag_branches ();
    void activate_edge (const uint32_t src_id, const uint32_t dest_id);
    void fill_branches (std::vector<double> &the_branches, const uint32_t num_branches);
    void compute_activation_times (const double cv);
    void compute_activation_times_with_diameter ();
    void compute_error (Graph *ref_network);
    void compute_error (std::vector<PMJ> pmjs);
    void compute_rmse_rrmse (Graph *ref_network, double &rmse, double &rrmse);
    void compute_rmse_rrmse (std::vector<PMJ> ref_pmjs, double &rmse, double &rrmse);
    void compute_min_max_lat (double &min_value, double &max_value);
    void compute_max_error (double &max_error, Graph *ref_network);
    void compute_max_error (double &max_error, std::vector<PMJ> ref_pmjs);
    void compute_epsilon_percentage (Graph *ref_network, const double epsilon, double &percentage);
    void compute_epsilon_percentage (std::vector<PMJ> ref_pmjs, const double epsilon, double &percentage);
    void compute_geometrics ();
    void convert_to_initial_network_topology (const char filename[]);
    void get_segment_length (std::vector<double> &the_segments);
    void get_bifurcation_angles (std::vector<double> &the_angles);
    void get_branches (std::vector<double> &the_branches);
    uint32_t get_root_position_index (const double root_pos[]);
    void print ();
    void write_terminals (const char filename[]);
    void write_active_pmjs (const char filename[], const double perc);
    void write_minimum_network (std::string pmj_filename, std::string filename);
    void write_graph_info (const char filename[]);
    void write_branches (const char filename[]);
    void write_LAT (const char filename[]);
    void write_polyline (const char filename[]);

    void build_and_sort_segments (std::vector< std::pair<uint32_t,uint32_t> > &sorted_segments);
    void build_shocker_structure_from_graph(std::vector< std::pair<uint32_t,uint32_t> > sorted_segments,\
                                            std::vector<Point_Segment*> &points, std::vector<Segment*> &segments);
    void build_shocker_points (std::vector<Point_Segment*> &points);
    void build_shocker_segments (std::vector< std::pair<uint32_t,uint32_t> > sorted_segments,\
                                std::vector<Point_Segment*> points, std::vector<Segment*> &segments);
    void update_diameter_and_LAT (std::vector<Segment*> &segments, const double diameter);
    void update_flags (std::vector<Point_Segment*> &points, std::vector<Segment*> &segments);
    void update_reference_closest_terminal_indexes (Graph *ref_network);
    void write_shocker_initial_network_topology(const char filename[], std::vector<Point_Segment*> points, std::vector<Segment*> segments);
};



#endif
