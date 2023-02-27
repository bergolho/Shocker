#ifndef UTILS_H
#define UTILS_H

#include <cstdio>
#include <cstdlib>
#include <cstdbool>
#include <cstdint>
#include <cstring>
#include <cmath>
#include <vector>
#include <algorithm>

#define UM_TO_MM 1.0e-03

class Point
{
public:
    Point (const uint32_t id, const double x, const double y, const double z)
    {
        this->id = id;
        this->x = x;
        this->y = y;
        this->z = z;
    }
public:
    uint32_t id;
    double x, y, z;
};

class Line
{
public:
    Line (const uint32_t src, const uint32_t dest)
    {
        this->src = src;
        this->dest = dest;
    }
public:
    uint32_t src;
    uint32_t dest;
};

class Cell
{
public:
    Cell (const uint32_t ids[])
    {
        for (uint32_t i = 0; i < 8; i++)
            this->ids[i] = ids[i];
    }
public:
    uint32_t ids[6];
};

void build_unitary_vector (double d[], const double src[], const double dest[]);
double calc_angle_between_vectors (const double u[], const double v[]);
double calc_norm (const double x1, const double y1, const double z1,\
                const double x2, const double y2, const double z2);
void compute_mean_std (std::vector<double> arr, double &mean, double &std);
void write_data_to_file (const char filename[], std::vector<double> arr);
void write_geometric_info_to_file (const char filename[], const double mean_segment_length, const double std_segment_length,\
                                    const double mean_branch_length, const double std_branch_length,\
                                    const double mean_bifurcation_angle, const double std_bifurcation_angle,\
                                    const uint32_t num_segments, const uint32_t num_branches, const uint32_t num_angles);
void write_electric_info_to_file(const char filename[], const double min_lat, const double max_lat, const double max_error,\
                                const double rmse, const double rrmse, const double eps2ms, const double eps5ms);
bool compareXYZ (Point p1, Point p2);

#endif