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
    Point (const uint32_t id, const double pos[])
    {
        this->id = id;
        memcpy(this->pos,pos,sizeof(double)*3);
    }
public:
    uint32_t id;
    double pos[3];
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


#endif