//
// Created by bergolho on 17/08/21.
//

#ifndef UTILS_H
#define UTILS_H

#include <cstdio>
#include <cstdlib>
#include <cstdbool>
#include <cstdint>
#include <cmath>
#include <sys/stat.h> 
#include <sys/types.h> 

#include <string>
#include <sstream>
#include <vector>
#include <numeric>
#include <algorithm>

// ==========================================================================================================================
// CONSTANTS AND MACROS
static const double COLLISION_THREASHOLD = 1.0E-04;     // Threashold size for considering that two segments intersect
static const std::string LINE_1 = "========================================================================";
static const std::string LINE_2 = "------------------------------------------------------------------------";
// ==========================================================================================================================

void usage (const char pname[]);
void usage_multiple_network (const char pname[]);
void create_directory (const char *path);
void print_message_and_exit (const char message[]);
void calc_mean_std (std::vector<double> arr, double &mean, double &std);
void write_vector_to_file (std::vector<double> arr, std::string filename);
void write_geometric_info_to_file (std::string filename,\
                    std::vector<double> segments, const double mean_segment_length, const double std_segment_length,\
                    std::vector<double> angles, const double mean_biff_angle, const double std_biff_angle);
void write_electric_info_to_file (std::string filename,\
                    const double max_lat_error, const double min_ref_lat, const double max_ref_lat,\
                    const double min_aprox_lat, const double max_aprox_lat,\
                    const double rmse, const double rrmse,\
                    const double epsilon_2ms, const double epsilon_5ms);

double euclidean_norm (const double x1, const double y1, const double z1,\
                    const double x2, const double y2, const double z2);
double calc_angle_between_vectors (const double u[], const double v[]);
double calc_dot_product (const double u[], const double v[]);
void calc_cross_product (const double a[], const double b[], double c[]);

double calc_lmin (const double l_d, const uint32_t num_terminals);

void calc_vector_subtract (const double a[], const double b[], double c[]);
double calc_vector_norm (const double a[]);
void normalize_vector (double arr[]);

//bool read_points_from_vtk (const char filename[], std::vector<Point*> &points);
//bool write_points_to_vtk (std::string filename, std::vector<Point*> input_points);
//bool read_pmjs_from_vtk (const char filename[], std::vector<PMJ*> &pmjs);
//bool write_pmjs_to_vtk (std::string filename, std::vector<PMJ*> pmjs);

bool check_size (const double p[]);
bool check_collinear (const double x1, const double y1, const double z1,\
                          const double x2, const double y2, const double z2,\
                          const double x3, const double y3, const double z3,\
                          const double x4, const double y4, const double z4);
bool collision_detection (const double x1, const double y1, const double z1,\
                          const double x2, const double y2, const double z2,\
                          const double x3, const double y3, const double z3,\
                          const double x4, const double y4, const double z4);
bool check_user_input (int argc, char *argv[]);
bool check_extension (std::string filename, std::string target_extension);
bool check_folder (std::string filename);

std::string get_color_code (const int color);
std::vector<uint32_t> argsort_abs (std::vector<double> arr);

#endif
