#ifndef GRAPH_H
#define GRAPH_H

#include <cstdio>
#include <cstdlib>
#include <cstdbool>
#include <cstdint>
#include <cstring>
#include <vector>
#include <queue>
#include <string>

#include "../utils/utils.h"
#include "../reader/reader.h"

#include "pqueue.h"

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
    uint32_t dest_id;
    double length;
public:
    Edge (const uint32_t dest_id, const double length)
    {
        this->dest_id = dest_id;
        this->length = length;
    }
};

class Node
{
public:
    uint32_t id;
    double pos[3];
    double sigma;
    std::vector<Edge> list_edges;
public:
    Node () { }
    Node (const uint32_t id, const double pos[], const double sigma)
    {
        this->id = id;
        memcpy(this->pos,pos,sizeof(double)*3);
        this->sigma = sigma;
    }
    void setNode (const uint32_t id, const double pos[], const double sigma)
    {
        this->id = id;
        memcpy(this->pos,pos,sizeof(double)*3);
        this->sigma = sigma;
    }
    void print ()
    {
        printf("|| %u (%g %g %g) [%g] ||",this->id,this->pos[0],this->pos[1],this->pos[2],this->sigma);
        for (uint32_t i = 0; i < this->list_edges.size(); i++)
            printf(" --> || %u %g || ",this->list_edges[i].dest_id,this->list_edges[i].length);
        printf("\n");
    }
};

class Graph
{
public:
    const double REF_CV = 1900.0;                   // Reference propagation velocity 
    uint32_t total_nodes;                           
    uint32_t total_edges;
    std::vector<Node> list_nodes;                   // List of nodes
    std::vector<double> dist;                       // Shortest distance of each cell from the root point
    std::vector<double> lat;                        // Local activation time of each cell
    std::vector<int> parent;                        // Parent indexes of each cell
    std::vector<uint32_t> terminals_indexes;        // Active PMJ's indexes
public:
    Graph ();
    Graph (std::vector<Point> points, std::vector<Line> lines);
    Graph (std::vector<Point> points, std::vector<Line> lines, std::vector<double> point_scalars);
    uint32_t get_closest_point (Point p);
    void depth_first_search (const uint32_t src_id, std::vector<double> &the_segments);
    void dfs (Node u, std::vector<bool> &dfs_visited, std::vector<double> &segments, double &segment_size, uint32_t &flag, uint32_t &total_segments);
    void dijkstra (const uint32_t src_id);
    void dijkstra_2 (const uint32_t src_id);
    void fill_terminal_indexes ();
    void compute_activation_times (const double cv);
    void compute_geometrics ();
    void get_segment_length (std::vector<double> &the_segments);
    void get_bifurcation_angles (std::vector<double> &the_angles);
    void print ();
    void write_terminals (const char filename[]);
    void write_network (const char filename[]);
    void write_LAT (const char filename[]);
    void write_MonoAlg3D (const char filename[]);
};

Graph* merge_networks (Graph *his, Graph *lv, Graph *rv, std::vector<Point> the_roots);

#endif
