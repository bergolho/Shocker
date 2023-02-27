#ifndef UTILS_H
#define UTILS_H

#include <cstdio>
#include <cstdlib>
#include <cstdbool>
#include <cstdint>
#include <cstring>
#include <cmath>
#include <vector>
#include <bitset>

class Point
{
public:
    Point (const uint32_t id, const double x, const double y, const double z)
    {
        this->id = id;
        this->x = x;
        this->y = y;
        this->z = z;
        this->value = 0.0;
    }
    Point (const uint32_t id, const double x, const double y, const double z, const double value)
    {
        this->id = id;
        this->x = x;
        this->y = y;
        this->z = z;
        this->value = value;
    }
    void print ()
    {
        printf("[Point %u] --> (%g %g %g), LAT = %g\n",this->id,this->x,this->y,this->z,this->value);
    }
public:
    uint32_t id;
    double x, y, z;
    double value;
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
    Cell (const uint32_t ids[], const double center[])
    {
        memcpy(this->ids,ids,sizeof(uint32_t)*8);
        memcpy(this->center,center,sizeof(double)*3);
    }
public:
    uint32_t ids[6];
    double center[3];
};

void build_unitary_vector (double d[], const double src[], const double dest[]);
double calc_angle_between_vectors (const double u[], const double v[]);
double calc_norm (const double x1, const double y1, const double z1,\
                const double x2, const double y2, const double z2);
void print_progress_bar (const uint32_t cur_iter, const uint32_t max_iter);                
void write_data_to_file (const char filename[], std::vector<double> arr);
void write_points_to_vtk (std::string filename, std::vector<Point> arr, const double scale);

#endif