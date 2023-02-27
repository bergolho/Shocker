#ifndef READER_H
#define READER_H

#include <cstdio>
#include <cstdlib>
#include <cstdbool>
#include <cstdint>
#include <cstring>
#include <cassert>
#include <string>
#include <unistd.h>
#include <ctype.h>

#include "../utils/utils.h"

class VTK_Reader
{
public:
    std::vector<Point> the_points;
    std::vector<Line> the_lines;
    std::vector<std::vector<double>> the_celldata;
public:
    VTK_Reader () { }
    VTK_Reader (std::string filename);
    uint32_t search_position (const double pos[]);
    void print ();
    void write (std::string filename);
};

class PMJ
{
public:
    uint32_t id;
    double x, y, z;
    double lat;
public:
    PMJ ();
    PMJ (const uint32_t id, const double x, const double y, const double z, const double lat = 0.0);
    void print ();
};

void read_points (std::string filename, std::vector<Point> &the_points);
void read_pmjs (std::string filename, std::vector<PMJ> &the_pmjs);
void read_root_position (std::string filename, double root_pos[]);
void compute_min_max_lat_pmjs (std::vector<PMJ> the_pmjs, double &min_value, double &max_value);
uint32_t get_closest_pmj (const double x, const double y, const double z, std::vector<PMJ> the_pmjs);

#endif