#include "reader.h"

VTK_Reader::VTK_Reader (std::string filename)
{
    char str[500];
    uint32_t np, nl, dummy;
    FILE *file = fopen(filename.c_str(),"r");
    if (!file)
    {
        fprintf(stderr,"[-] ERROR! Cannot read file '%s' !\n",filename.c_str());
        exit(1);
    }

    while (fscanf(file,"%s",str) != EOF)
        if (strcmp(str,"POINTS") == 0) break;
    fscanf(file,"%u %s",&np,str);
    for (uint32_t i = 0; i < np; i++)
    {
        double pos[3];
        fscanf(file,"%lf %lf %lf",&pos[0],&pos[1],&pos[2]);

        Point p(i,pos[0],pos[1],pos[2]);
        this->the_points.push_back(p);
    }
    while (fscanf(file,"%s",str) != EOF)
        if (strcmp(str,"LINES") == 0) break;
    fscanf(file,"%u %u",&nl,&dummy);
    for (uint32_t i = 0; i < nl; i++)
    {
        uint32_t ids[2];
        fscanf(file,"%u %u %u",&dummy,&ids[0],&ids[1]);

        if (ids[0] != ids[1])
        {
            Line l(ids[0],ids[1]);
            this->the_lines.push_back(l);
        }
    }

    // Read the CELL_DATA section, if it exist
    if (fscanf(file,"%s",str) != EOF) {
        
        if (strcmp(str,"CELL_DATA") == 0) {

            uint32_t num_fields;
            while (fscanf(file,"%s",str) != EOF) {
                if (strcmp(str,"FieldData") == 0) break;
            }
            fscanf(file,"%u",&num_fields);
            the_celldata.assign(num_fields,std::vector<double>());

            for (uint32_t k = 0; k < num_fields; k++) {

                while (fscanf(file,"%s",str) != EOF) {
                    if (strcmp(str,"float") == 0) break;
                }

                for (int i = 0; i < nl; i++) {
                    double scalar;
                    fscanf(file,"%lf",&scalar);
                    the_celldata[k].push_back(scalar);
                }
            }
        }
    }

    //print();
}

uint32_t VTK_Reader::search_position (const double pos[])
{
    uint32_t closest_index = 0;
    double closest_dist = __DBL_MAX__;
    for (uint32_t i = 0; i < this->the_points.size(); i++)
    {
        double x = this->the_points[i].x;
        double y = this->the_points[i].y;
        double z = this->the_points[i].z;

        double dist = sqrt( pow(x-pos[0],2) + pow(y-pos[1],2) + pow(z-pos[2],2) );
        if (dist < closest_dist)
        {
            closest_dist = dist;
            closest_index = i;
        }
    }
    return closest_index;
}

void read_points (std::string filename, std::vector<Point> &the_points)
{
    uint32_t np;
    FILE *file = fopen(filename.c_str(),"r");
    fscanf(file,"%u",&np);
    for (uint32_t i = 0 ; i < np; i++)
    {
        double x, y, z;
        fscanf(file,"%lf %lf %lf",&x,&y,&z);
        Point p(i,x,y,z);
        the_points.push_back(p);
    }
    fclose(file);
}

void read_pmjs_from_vtk (std::string filename, std::vector<Point> &the_pmjs)
{

}

void VTK_Reader::print ()
{
    printf("======= POINTS =======\n");
    for (uint32_t i = 0; i < this->the_points.size(); i++)
    {
        uint32_t id = this->the_points[i].id;
        double x = this->the_points[i].x;
        double y = this->the_points[i].y;
        double z = this->the_points[i].z;

        printf("[%u] = (%g %g %g)\n",id,x,y,z);
    }

    printf("======= LINES =======\n");
    for (uint32_t i = 0; i < this->the_lines.size(); i++)
    {
        uint32_t id = i;
        uint32_t src = this->the_lines[i].src;
        uint32_t dest = this->the_lines[i].dest;

        printf("[%u] %u --> %u\n",i,src,dest);
    }

    if (this->the_celldata.size() != 0) {
        printf("======= CELL_DATA =======\n");
        for (uint32_t k = 0; k < this->the_celldata.size(); k++) {
            printf("CELL_DATA field %u\n",k);
            for (uint32_t i = 0; i < this->the_celldata[k].size(); i++) {
                printf("[%u] = %g\n",i,this->the_celldata[k][i]);
            }
        }
        printf("\n");
    }
}

void VTK_Reader::write (std::string filename)
{
    FILE *file = fopen(filename.c_str(),"w+");
    fprintf(file,"# vtk DataFile Version 4.1\n");
    fprintf(file,"vtk output\n");
    fprintf(file,"ASCII\n");
    fprintf(file,"DATASET POLYDATA\n");
    fprintf(file,"POINTS %u float\n",this->the_points.size());
    for (uint32_t i = 0; i < this->the_points.size(); i++)
        fprintf(file,"%g %g %g\n",this->the_points[i].x,this->the_points[i].y,this->the_points[i].z);
    fprintf(file,"LINES %u %u\n",this->the_lines.size(),this->the_lines.size()*3);
    for (uint32_t i = 0; i < this->the_lines.size(); i++)
        fprintf(file,"2 %u %u\n",this->the_lines[i].src,this->the_lines[i].dest);
    fclose(file);
}

PMJ::PMJ ()
{
    this->id = 0;
    this->x = 0.0;
    this->y = 0.0;
    this->z = 0.0;
    this->lat = 0.0;
}

PMJ::PMJ (const uint32_t id, const double x, const double y, const double z, const double lat)
{
    this->id = id;
    this->x = x;
    this->y = y;
    this->z = z;
    this->lat = lat;
}

void PMJ::print ()
{
    printf("PMJ %u = ( %g %g %g ) || LAT = %g\n",this->id,this->x,this->y,this->z,this->lat);
}

void read_pmjs (std::string filename, std::vector<PMJ> &the_pmjs)
{
    uint32_t np;
    char str[200];
    FILE *file = fopen(filename.c_str(),"r");
    while (fscanf(file,"%s",str) != EOF)
        if (strcmp(str,"POINTS") == 0) break;
    fscanf(file,"%u %s",&np,str);
    for (uint32_t i = 0; i < np; i++)
    {
        double x, y, z;
        fscanf(file,"%lf %lf %lf",&x,&y,&z);
        PMJ pmj(i,x,y,z);
        the_pmjs.push_back(pmj);
    }
    while (fscanf(file,"%s",str) != EOF)
        if (strcmp(str,"float") == 0) break;
    for (uint32_t i = 0; i < np; i++)
    {
        double lat;
        fscanf(file,"%lf",&lat);
        the_pmjs[i].lat = lat;
    }
    fclose(file);
}

void compute_min_max_lat_pmjs (std::vector<PMJ> the_pmjs, double &min_value, double &max_value)
{
    uint32_t num_pmjs = the_pmjs.size();
    max_value = __DBL_MIN__;
    min_value = __DBL_MAX__;
    for (uint32_t i = 0; i < num_pmjs; i++)
    {
        if (the_pmjs[i].lat > max_value)
        {
            max_value = the_pmjs[i].lat;
        }
        if (the_pmjs[i].lat < min_value)
        {
            min_value = the_pmjs[i].lat;
        }
    }
}

uint32_t get_closest_pmj (const double x, const double y, const double z, std::vector<PMJ> the_pmjs)
{
    uint32_t closest_id = 0;
    double closest_dist = __DBL_MAX__;
    for (uint32_t j = 0; j < the_pmjs.size(); j++)
    {
        double dist = calc_norm(the_pmjs[j].x,the_pmjs[j].y,the_pmjs[j].z,x,y,z);
        if (dist < closest_dist)
        {
            closest_dist = dist;
            closest_id = j;
        }
    }
    return closest_id;
}

void read_root_position (std::string filename, double root_pos[])
{
    FILE *file = fopen(filename.c_str(),"r");
    if (!file)
    {
        fprintf(stderr,"[-] ERROR! Cannot read filename '%s'!\n",filename.c_str());
        exit(EXIT_FAILURE);
    }
    fscanf(file,"%lf %lf %lf",&root_pos[0],&root_pos[1],&root_pos[2]);
    fclose(file);
}