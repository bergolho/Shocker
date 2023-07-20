//
// Created by bergolho on 18/08/21.
//

#ifndef CLOUD_H
#define CLOUD_H

#include <cstdio>
#include <cstdlib>
#include <cstdbool>
#include <cstdint>
#include <cmath>
#include <cassert>

#include <vector>
#include <set>
#include <map>
#include <string>
#include <algorithm>

#include <vtkXMLPolyDataReader.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkPolyDataMapper.h>
#include <vtkGenericDataObjectReader.h>
#include <vtkProperty.h>
#include <vtkPolyData.h>
#include <vtkUnstructuredGrid.h>
#include <vtkIdList.h>
#include <vtkSmartPointer.h>
#include <vtkPolyDataMapper.h>
#include <vtkDijkstraGraphGeodesicPath.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkLine.h>
#include <vtkPointLocator.h>
#include <vtkTriangle.h>

#include "../../options/cloud_config.h"
#include "../../options/pmj_config.h"
#include "../../utils/utils.h"

#include "point.h"
#include "pmj.h"

class Cloud_Data
{
public:
    static constexpr double MIN_REGION_RADIUS = 10;         // {um}
    static constexpr double REGION_RADIUS_OFFSET = 50;      // {um}
    static constexpr double REGION_RADIUS = 7500;            // {um}
    uint32_t cur_id;
    uint32_t rand_offset;
    double phi;
    std::string cloud_filename;
    std::vector<CloudPoint*> points;
    uint32_t counter_cloud_passes = 0;
public:
    Cloud_Data (CloudConfig *cloud_config, std::string output_dir, const double root_pos[]);
    ~Cloud_Data ();
    CloudPoint* sort_point ();
    void pre_process (const double root_pos[]);
    void calc_distal_root_position (const double x_prox[], const double l_d, double x_dist[]);
    void print ();
    void write (std::string output_dir);
private:
    bool read_cloud_from_vtk ();
    void remap (const double root_pos[]);
    void filter ();

};

class PMJ_Data
{
public:
    uint32_t total_num_connected;                           // Total number of connected PMJ's within the LAT tolerance
    uint32_t connection_rate;                               // Connection rate to process PMJs
    double lat_error_tolerance;                             // LAT error tolerance
    std::string pmj_filename;                               // Filename with the PMJs location
    std::vector<PMJPoint*> pmjs;                            // Array with the PMJs
public:
    PMJ_Data (PMJConfig *config, std::string output_dir);
    ~PMJ_Data ();
    bool has_pmjs ();
    void reset ();
    void print ();
    void write (std::string output_dir);
    void write_errors (std::string filename);
    void write_errors_in_vtk (std::string filename);
private:
    bool read_pmjs_from_vtk ();
    void filter ();
};

class Surface_Data
{
public:
    vtkSmartPointer<vtkPolyData> endo_data;               
    vtkSmartPointer<vtkDijkstraGraphGeodesicPath> dijkstra;     
    vtkSmartPointer<vtkPointLocator> point_locator;             
public:
    Surface_Data (std::string filename);
    ~Surface_Data ();
    void read_surface_from_file (std::string filename);
    std::vector<CloudPoint> compute_geodesic_pathway (const double a[], const double b[], const bool use_first_node, const double m[]);
    void print ();
};

class Cloud
{
public:
    Cloud_Data *cloud_data;
    PMJ_Data *pmj_data;
    Surface_Data *surface_data;
public:
    Cloud (CloudConfig *cloud_config, PMJConfig *pmj_config, std::string output_dir, const double root_pos[]);
    ~Cloud ();
    void tag_pmjs_to_cloud ();
    void print ();
    void write (std::string output_dir);
};

void write_cloud_array_to_file (std::string filename, std::vector<CloudPoint*> arr);

#endif