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