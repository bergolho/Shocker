#ifndef READER_H
#define READER_H

#include <cstdio>
#include <cstdlib>
#include <cstdbool>
#include <cstdint>
#include <cstring>
#include <unistd.h>
#include <ctype.h>

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkUnstructuredGrid.h>
#include <vtkCell.h>
#include <vtkLine.h>
#include <vtkHexahedron.h>
#include <vtkFloatArray.h>
#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkGenericDataObjectReader.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkPolyDataWriter.h>

#include "../pmj/pmj.h"
#include "../utils/utils.h"

#define REGION_RADIUS 5000.0

class VTU_Reader
{
public:
    std::vector<Point> the_points;
    std::vector<Cell> the_cells;
    std::vector<double> the_scalars;
public:
    VTU_Reader () { }
    VTU_Reader (std::string filename);
    void get_cells_center_positions(std::vector<Point> &out);
    void print ();
};

void get_closest_points (std::vector<Point> in, std::vector<Point> in2, std::vector<Point> &out, const double percentage, const uint32_t offset);
void read_points_from_vtk (std::string filename, std::vector<Point> &arr);
void read_pmjs (std::string filename, std::vector<PMJ> &pmjs);

std::vector<Point> filter_points (std::vector<Point> endo, const uint32_t target_num_points);
std::vector<Point> filter_points_with_file (std::vector<Point> endo, std::vector<Point> pmjs, const uint32_t target_num_points);
void concatenate_points (std::vector<Point> &out, std::vector<Point> in);

#endif
