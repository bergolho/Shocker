#include "utils.h"

void create_directory (const char *path)
{
    if (mkdir(path,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) != -1)
        printf("[INFO] Output directory created at:> %s\n",path);
    else
        fprintf(stderr,"[INFO] Output directory already exists!\n");
}

void print_message_and_exit (const char message[])
{
    printf("[INFO] Message:> %s\n",message);
    exit(EXIT_SUCCESS);
}

void write_vector_to_file (std::vector<double> arr, std::string filename)
{
    FILE *file = fopen(filename.c_str(),"w+");
    for (uint32_t i = 0; i < arr.size(); i++) fprintf(file,"%g\n",arr[i]);
    fclose(file);
}

void write_geometric_info_to_file (std::string filename,\
                    std::vector<double> segments, const double mean_segment_length, const double std_segment_length,\
                    std::vector<double> angles, const double mean_biff_angle, const double std_biff_angle)
{
    FILE *file = fopen(filename.c_str(),"w+");
    fprintf(file,"[INFO] Total number of segment = %lu\n",segments.size());
    fprintf(file,"[INFO] Segment length = %g +/- %g mm\n",mean_segment_length,std_segment_length);
    fprintf(file,"[INFO] Total number of bifurcations = %lu\n",angles.size());
    fprintf(file,"[INFO] Bifurcation angle = %g +/- %g degrees\n",mean_biff_angle,std_biff_angle);
    fprintf(file,"%.2lf +/- %.2lf\t%.2lf +/- %.2lf\t%lu\t%lu\n",mean_segment_length,std_segment_length,\
                                                                mean_biff_angle,std_biff_angle,
                                                                segments.size(),angles.size());
    fclose(file);
}

void write_electric_info_to_file (std::string filename,\
                    const double max_lat_error, const double min_ref_lat, const double max_ref_lat,\
                    const double min_aprox_lat, const double max_aprox_lat,\
                    const double rmse, const double rrmse,\
                    const double epsilon_2ms, const double epsilon_5ms)
{
    FILE *file = fopen(filename.c_str(),"a");
    fprintf(file,"[INFO] Reference --> || min.LAT = %.2lf ms || max.LAT = %.2lf ms ||\n",min_ref_lat,max_ref_lat);
    fprintf(file,"[INFO] Aproximation --> || min.LAT = %.2lf ms || max.LAT = %.2lf ms || max.ERROR = %.2lf ms\n",min_aprox_lat,max_aprox_lat,max_lat_error);
    fprintf(file,"[INFO] || RMSE = %.2lf ms || RRMSE = %.2lf %% ||\n",rmse,rrmse);
    fprintf(file,"[INFO] || Epsilon 2ms = %.2lf %% || Epsilon 5ms = %.2lf ||\n",epsilon_2ms,epsilon_5ms);
    fprintf(file,"%.2lf\t%.2lf\t%.2lf\t%.2lf\t%.2lf\t%.2lf\t%.2lf\n",min_aprox_lat,max_aprox_lat,max_lat_error,rmse,rrmse,epsilon_2ms,epsilon_5ms);
    fclose(file);
}

void calc_mean_std (std::vector<double> arr, double &mean, double &std)
{
    double sum = 0.0;
    double n = (double)arr.size();

    for (uint32_t i = 0; i < arr.size(); i++) sum += arr[i];
    mean = sum / n;

    sum = 0.0;
    for (uint32_t i = 0; i < arr.size(); i++) sum += pow(arr[i]-mean,2);
    std = sqrt(sum / n);
}

void normalize_vector (double arr[])
{
    double norm = sqrt(pow(arr[0],2)+pow(arr[1],2)+pow(arr[2],2));
    for (int i = 0; i < 3; i++)
        arr[i] /= norm;
}

void calc_vector_subtract (const double a[], const double b[], double c[])
{
    for (int i = 0; i < 3; i++)
        c[i] = a[i] - b[i];
}

double calc_vector_norm (const double a[])
{
    return sqrt(pow(a[0],2)+pow(a[1],2)+pow(a[2],2));
}

std::vector<uint32_t> argsort_abs (std::vector<double> arr)
{
    std::vector<uint32_t> sorted_indexes(arr.size());
    std::iota(sorted_indexes.begin(),sorted_indexes.end(),0);
    std::sort(sorted_indexes.begin(),sorted_indexes.end(), [&](int i,int j) { return fabs(arr[i]) < fabs(arr[j]); });
    return sorted_indexes;
}

double euclidean_norm (const double x1, const double y1, const double z1,\
                    const double x2, const double y2, const double z2)
{
    return sqrt( pow(x2-x1,2) + pow(y2-y1,2) + pow(z2-z1,2) );
}

double calc_angle_between_vectors (const double u[], const double v[])
{
    double dot_product = u[0]*v[0] + u[1]*v[1] + u[2]*v[2];

    double angle_radians = acos(dot_product);
    if (std::isnan(angle_radians)) angle_radians = 0.0;

    // Return the angle in degrees
    return angle_radians * 180.0 / M_PI;
}

double calc_dot_product (const double u[], const double v[])
{
    return (u[0] * v[0]) + (u[1] * v[1]) + (u[2] * v[2]); 
}

void calc_cross_product (const double a[], const double b[], double c[])
{
    c[0] = a[1] * b[2] - a[2] * b[1];
    c[1] = a[2] * b[0] - a[0] * b[2];
    c[2] = a[0] * b[1] - a[1] * b[0];
}

double calc_lmin (const double l_d, const uint32_t num_terminals)
{
    return sqrt( (l_d*l_d)/num_terminals );
}

bool check_size (const double p[])
{
    return (fabs(p[0]) < COLLISION_THREASHOLD && fabs(p[1]) < COLLISION_THREASHOLD && fabs(p[2]) < COLLISION_THREASHOLD) ? true : false;
}

bool check_collinear (const double x1, const double y1, const double z1,\
                          const double x2, const double y2, const double z2,\
                          const double x3, const double y3, const double z3,\
                          const double x4, const double y4, const double z4)
{
    if (x1 == x3 && y1 == y3 && z1 == z3 &&\
        x2 == x4 && y2 == y4 && z2 == z4)
        return true;
    else
        return false;
}

// Return true if we detect a collision
/*
   Calculate the line segment PaPb that is the shortest route between
   two lines P1P2 and P3P4. Calculate also the values of mua and mub where
      Pa = P1 + mua (P2 - P1)
      Pb = P3 + mub (P4 - P3)
   Return FALSE if no solution exists.
*/
bool collision_detection (const double x1, const double y1, const double z1,\
                          const double x2, const double y2, const double z2,\
                          const double x3, const double y3, const double z3,\
                          const double x4, const double y4, const double z4)
{
    double p13[3], p43[3], p21[3];

    if (check_collinear(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4)) return true;

    p13[0] = x1 - x3;
    p13[1] = y1 - y3;
    p13[2] = z1 - z3;

    p43[0] = x4 - x3;
    p43[1] = y4 - y3;
    p43[2] = z4 - z3;

    p21[0] = x2 - x1;
    p21[1] = y2 - y1;
    p21[2] = z2 - z1;
    
    if (check_size(p43)) return false;
    if (check_size(p21)) return false;

    double dot_product_1343 = calc_dot_product(p13,p43);
    double dot_product_4321 = calc_dot_product(p43,p21);
    double dot_product_1321 = calc_dot_product(p13,p21);
    double dot_product_4343 = calc_dot_product(p43,p43);
    double dot_product_2121 = calc_dot_product(p21,p21);

    double numerator = dot_product_1343 * dot_product_4321 - dot_product_1321 * dot_product_4343;
    double denominator = dot_product_2121 * dot_product_4343 - dot_product_4321 * dot_product_4321;
    if (fabs(denominator) < COLLISION_THREASHOLD) return false;
    
    double mua = numerator / denominator;
    double mub = (dot_product_1343 + dot_product_4321 * mua) / dot_product_4343;

    double pa[3], pb[3];
    pa[0] = x1 + mua * p21[0];
    pa[1] = y1 + mua * p21[1];
    pa[2] = z1 + mua * p21[2];

    pb[0] = x3 + mub * p43[0];
    pb[1] = y3 + mub * p43[1];
    pb[2] = z3 + mub * p43[2];

    if ( (mua > 0.0 && mua < 1) && (mub > 0.0 && mub < 1.0) )
    {
        double norm = sqrt( pow(pa[0]-pb[0],2) + pow(pa[1]-pb[1],2) + pow(pa[2]-pb[2],2) );
        if (norm > COLLISION_THREASHOLD)
            return false;
        else
            return true;
    }
    else
        return false;
}

bool check_user_input (int argc, char *argv[])
{
    int num_networks = argc-1;
    // Check if all the input arguments are INI configuration files
    for (int i = 0; i < argc-2; i++)
    {
        std::string config_filename = argv[i+1];
        printf("%s\n",config_filename.c_str());
        if (!check_extension(config_filename,"ini")) return false;
    }
    // Check if the last argument is a directory
    return (check_folder(argv[argc-1])) ? true : false;
}

bool check_extension (std::string filename, std::string target_extension)
{
    size_t nlen = filename.size();
    std::string extension = filename.substr( nlen-3, nlen );
    return (extension == target_extension) ? true : false;
}

bool check_folder (std::string filename)
{
    size_t nlen = filename.size();
    size_t npos = filename.find_first_of(".");
    return (npos != nlen-4) ? true : false;
}

std::string get_color_code (const int color)
{
    switch (color)
    {
        case 0: return "\033[22;37m";
        case 1: return "\033[1;32m";
        case 2: return "\033[1;31m";
        case 3: return "\033[1;33m";
        case 4: return "\033[1;35m";
        case 5: return "\033[0;36m";
        default: return "\033[22;37m";
    }
}

void usage (const char pname[])
{
    printf("%s\n", LINE_2.c_str());
    printf("Usage:> %s <input_config>\n",pname);
    printf("%s\n", LINE_2.c_str());
    printf("Example:> %s inputs/simplified_LV.ini\n",pname);
    printf("%s\n", LINE_2.c_str());        
}

// TODO: Future projects ...
void usage_multiple_network (const char pname[])
{
    printf("%s\n", LINE_2.c_str());
    printf("Usage:> %s [list_of_input_config] <output_dir>\n",pname);
    printf("%s\n", LINE_2.c_str());
    printf("Example:> %s inputs/benchmark_3d_cuboid_north.ini inputs/benchmark_3d_cuboid_east.ini inputs/benchmark_3d_cuboid_west.ini outputs/benchmark3d:cuboid\n",pname);
    printf("%s\n", LINE_2.c_str());
}