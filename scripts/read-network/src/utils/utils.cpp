#include "utils.h"

double calc_norm (const double x1, const double y1, const double z1,\
                const double x2, const double y2, const double z2)
{
    return sqrt( pow(x1-x2,2) + pow(y1-y2,2) + pow(z1-z2,2) );
}

void build_unitary_vector (double d[], const double src[], const double dest[])
{
    double norm = calc_norm(src[0],src[1],src[2],\
                        dest[0],dest[1],dest[2]);
    if (norm < 1.0e-05)
        norm = 1.0e-05;

    for (uint32_t i = 0; i < 3; i++)
    {
        d[i] = (dest[i] - src[i]) / norm;
    }
        
}

double calc_angle_between_vectors (const double u[], const double v[])
{
    double dot_product = u[0]*v[0] + u[1]*v[1] + u[2]*v[2];

    if (dot_product > 1.0)
        dot_product = 1.0;
    else if (dot_product < -1.0)
        dot_product = -1.0;

    double angle_radians = acos(dot_product);
    if (dot_product == 1.0)
        angle_radians = 0.0;
    else if (dot_product == -1.0)
        angle_radians = M_PI;

    // Return the angle in degrees
    return angle_radians * 180.0 / M_PI;
}

void write_data_to_file (const char filename[], std::vector<double> arr)
{
    FILE *file = fopen(filename,"w+");

    for (uint32_t i = 0; i < arr.size(); i++)
        fprintf(file,"%g\n",arr[i]);

    fclose(file);
}

void compute_mean_std (std::vector<double> arr, double &mean, double &std)
{
    mean = 0.0;
    for (uint32_t i = 0 ; i < arr.size(); i++)
        mean += arr[i];
    mean /= (double)arr.size();

    std = 0.0;
    for (uint32_t i = 0 ; i < arr.size(); i++)
        std += powf(arr[i]-mean,2);
    std /= (double)arr.size();
    std = sqrt(std);
}

bool compareXYZ (Point p1, Point p2)
{
    return (p1.x < p2.x) || \
               ((!(p2.x < p1.x)) && (p1.y < p2.y)) || \
               ((!(p2.x < p1.x)) && (!(p2.y < p1.y)) && (p1.z < p2.z));
}

void write_geometric_info_to_file (const char filename[], const double mean_segment_length, const double std_segment_length,\
                                    const double mean_branch_length, const double std_branch_length,\
                                    const double mean_bifurcation_angle, const double std_bifurcation_angle,\
                                    const uint32_t num_segments, const uint32_t num_branches, const uint32_t num_angles)
{
    FILE *file = fopen(filename,"w+");

    // [Vertical output]
    // mean_segment_length, std_segment_length, mean_branch_length, std_branch_length, mean_biff_angle, std_biff_angle, num_segments, num_angles
    //fprintf(file,"%.2lf +/- %.2lf\n%.2lf +/- %.2lf\n%.2lf +/- %.2lf\n%u\n%u\n%u\n",mean_segment_length*UM_TO_MM,std_segment_length*UM_TO_MM,\
                                            mean_branch_length*UM_TO_MM,std_branch_length*UM_TO_MM,\
                                            mean_bifurcation_angle,std_bifurcation_angle,\
                                            num_segments,num_branches,num_angles);
    // [Horizontal output]
    // mean_segment_length, std_segment_length, mean_branch_length, std_branch_length, mean_biff_angle, std_biff_angle, num_segments, num_angles
    fprintf(file,"%.2lf +/- %.2lf\t%.2lf +/- %.2lf\t%.2lf +/- %.2lf\t%u\t%u\t%u\n",mean_segment_length*UM_TO_MM,std_segment_length*UM_TO_MM,\
                                            mean_branch_length*UM_TO_MM,std_branch_length*UM_TO_MM,\
                                            mean_bifurcation_angle,std_bifurcation_angle,\
                                            num_segments,num_branches,num_angles);

    fclose(file);
}

void write_electric_info_to_file(const char filename[], const double min_lat, const double max_lat, const double max_error,\
                                const double rmse, const double rrmse, const double eps2ms, const double eps5ms)
{
    FILE *file = fopen(filename,"w+");

    // minLAT, maxLAT, maxERROR, RMSE, RRMSE, epsilon2ms, epsilon5ms
    fprintf(file,"%.2lf\t%.2lf\t%.2lf\t%.2lf\t%.2lf\t%.2lf\t%.2lf\n",min_lat,max_lat,max_error,rmse,rrmse*100,eps2ms,eps5ms);

    fclose(file);
}