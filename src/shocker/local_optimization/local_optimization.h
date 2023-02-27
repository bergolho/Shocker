#ifndef LOCAL_OPTIMIZATION_H
#define LOCAL_OPTIMIZATION_H

#include <cstdio>
#include <cstdlib>
#include <cmath>

#include <string>
#include <vector>

#include "../../options/local_optimization_config.h"
#include "../network/segment.h"

class BifurcationPoint
{
public:
    uint32_t id;
    double pos[3];
public:
    BifurcationPoint (const uint32_t id, const double pos[]);
    ~BifurcationPoint ();
    inline void setId (const uint32_t id) { this->id = id; }
    inline void setPosition (const double pos[]) { memcpy(this->pos,pos,sizeof(double)*3); }
    void print ();
};

class LocalOptimization
{
public:
    static constexpr double NE = 10;             // Discretization of the bifurcation plane (Default = 5)
    std::vector<BifurcationPoint> biff_points;  // Array with the points inside the bifurcation plane
public:
    LocalOptimization (LocalOptimizationConfig *local_opt_config);
    ~LocalOptimization ();
    void update_bifurcation_points (Segment *iconn, Segment *ibiff, Segment *inew);
    void print ();
    void write (std::string output_dir, const double translate[]);
};

#endif