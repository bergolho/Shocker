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
    if (norm < 1.0e-08)
        norm = 1.0e-08;

    for (uint32_t i = 0; i < 3; i++)
    {
        d[i] = (dest[i] - src[i]) / norm;
    }
        
}

double calc_angle_between_vectors (const double u[], const double v[])
{
    double dot_product = u[0]*v[0] + u[1]*v[1] + u[2]*v[2];

    double angle_radians = acos(dot_product);

    // Return the angle in degrees
    return angle_radians * 180.0 / M_PI;
}

void print_progress_bar (const uint32_t cur_iter, const uint32_t max_iter)
{
    double percentage = (cur_iter / (double) max_iter) * 100.0;
    uint32_t filled_length = nearbyint(100 * cur_iter / max_iter);

    std::string the_bar;
    for (uint32_t i = 0; i < filled_length; i++)
    the_bar += "\u2588";    // Using Unicode here ... (filled box)
    for (uint32_t i = 0; i < 100-filled_length; i++)
    the_bar += "-";

    printf("%s |%s| %.1lf%% %s\r","Progress",the_bar.c_str(),percentage,"Complete"); // Carriage return
    fflush(stdout); // Flush the standard output
}

void write_data_to_file (const char filename[], std::vector<double> arr)
{
    FILE *file = fopen(filename,"w+");

    for (uint32_t i = 0; i < arr.size(); i++)
        fprintf(file,"%g\n",arr[i]);

    fclose(file);
}

void write_points_to_vtk (std::string filename, std::vector<Point> arr, const double scale)
{
    uint32_t np = arr.size();
    FILE *file = fopen(filename.c_str(),"w+");
    fprintf(file,"# vtk DataFile Version 4.1\n");
    fprintf(file,"vtk output\n");
    fprintf(file,"ASCII\n");
    fprintf(file,"DATASET POLYDATA\n");
    fprintf(file,"POINTS %u float\n",np);
    for (uint32_t i = 0; i < np; i++)
        fprintf(file,"%g %g %g\n",arr[i].x*scale,arr[i].y*scale,arr[i].z*scale);
    fprintf(file,"VERTICES %u %u\n",np,np*2);
    for (uint32_t i = 0; i < np; i++)
        fprintf(file,"1 %u\n",i);
    fprintf(file,"POINT_DATA %u\n",np);
    fprintf(file,"FIELD FieldData 1\n");
    fprintf(file,"LAT 1 %u float\n",np);
    for (uint32_t i = 0; i < np; i++)
        fprintf(file,"%g\n",arr[i].value);
    fclose(file);
}