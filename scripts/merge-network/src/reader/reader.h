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
    bool has_point_data;
    std::vector<Point> the_points;
    std::vector<Line> the_lines;
    std::vector<double> the_point_scalars;
public:
    VTK_Reader () { }
    VTK_Reader (std::string filename);
    uint32_t search_position (const double target_pos[]);
    void print ();
    void write (std::string filename);
};

void read_root_positions(std::string filename, std::vector<Point> &the_roots);
void read_points (std::string filename, std::vector<Point> &the_points);

#endif