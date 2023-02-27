#ifndef PMJ_H
#define PMJ_H

#include <cstdio>
#include <cstdlib>
#include <cstdbool>
#include <cstdint>
#include <cstring>

#include "../utils/utils.h"

#define PMJ_REGION_RADIUS 5000

class PMJ
{
public:
    bool is_active;             // Flag for activation
    uint32_t id;                // Index of the cell in the mesh
    uint32_t term_id;           // Index of the terminal
    double pos[3];              // Cordinates {um,um,um}
    double pk_lat;              // LAT of the Purkinje cell {ms}
    double tiss_lat;            // Mean LAT of the tissue cells {ms}
    double delay;               // PMJ delay {ms}
    double pk_dist;             // Distance from the root to the terminal {um}
    double pk_cv;               // Conduction velocity {um/ms}
public:
    PMJ ();
    PMJ (const uint32_t id, const uint32_t term_id, const double pos[],\
         const double pk_lat, const double tiss_lat, const double delay,\
         const double pk_dist, const double pk_cv);
    void print ();
};

void get_closest_points_inside_sphere (std::vector<Point> endo_points, std::vector<PMJ> pmj_points, std::vector<Point> &out, const double percentage, const double his_offset);

#endif
