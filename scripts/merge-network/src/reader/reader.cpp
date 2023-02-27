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

        Point p(i,pos);
        this->the_points.push_back(p);
    }
    while (fscanf(file,"%s",str) != EOF)
        if (strcmp(str,"LINES") == 0) break;
    fscanf(file,"%u %u",&nl,&dummy);
    for (uint32_t i = 0; i < nl; i++)
    {
        uint32_t ids[2];
        fscanf(file,"%u %u %u",&dummy,&ids[0],&ids[1]);

        Line l(ids[0],ids[1]);
        this->the_lines.push_back(l);
    }
    // Check for the POINT_DATA sectiom
    if (fscanf(file,"%s",str) != EOF) 
    {
        if (strcmp(str,"POINT_DATA") == 0) 
        {
            while (fscanf(file,"%s",str) != EOF)
                if (strcmp(str,"default") == 0) break;

            for (uint32_t i = 0; i < np; i++) 
            {
                double ps;
                fscanf(file,"%lf",&ps);
                this->the_point_scalars.push_back(ps);
            }
            this->has_point_data = true;
        }
    }
    else
    {
        this->has_point_data = false;
    }
}

uint32_t VTK_Reader::search_position (const double target_pos[])
{
    uint32_t closest_index = 0;
    double closest_dist = __DBL_MAX__;
    for (uint32_t i = 0; i < this->the_points.size(); i++)
    {
        double *pos = this->the_points[i].pos;
        double dist = sqrt( pow(pos[0]-target_pos[0],2) + pow(pos[1]-target_pos[1],2) + pow(pos[2]-target_pos[2],2) );
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
        double pos[3];
        fscanf(file,"%lf %lf %lf",&pos[0],&pos[1],&pos[2]);
        Point p(i,pos);
        the_points.push_back(p);
    }
    fclose(file);
}

void read_root_positions(std::string filename, std::vector<Point> &the_roots)
{
    double pos[3];
    FILE *file = fopen(filename.c_str(),"r");
    if (!file)
    {
        fprintf(stderr,"[-] ERROR! Cannot read root position file '%s'\n",filename.c_str());
        exit(EXIT_FAILURE);
    }
    for (uint32_t i = 0; i < 2; i++)
    {
        if (!fscanf(file,"%lf %lf %lf",&pos[0],&pos[1],&pos[2]))
        {
            fprintf(stderr,"[-] ERROR! Reading root positions file!\n");
            exit(1);
        }
        Point p(i,pos);
        the_roots.push_back(p);
    }
    fclose(file);

    assert( the_roots.size() == 2 );
}

void VTK_Reader::print ()
{
    printf("======= POINTS =======\n");
    for (uint32_t i = 0; i < this->the_points.size(); i++)
    {
        uint32_t id = this->the_points[i].id;
        double *pos = this->the_points[i].pos;

        printf("[%u] = (%g %g %g)\n",id,pos[0],pos[1],pos[2]);
    }
    printf("======= LINES =======\n");
    for (uint32_t i = 0; i < this->the_lines.size(); i++)
    {
        uint32_t id = i;
        uint32_t src = this->the_lines[i].src;
        uint32_t dest = this->the_lines[i].dest;

        printf("[%u] %u --> %u\n",i,src,dest);
    }
}

void VTK_Reader::write (std::string filename)
{
    uint32_t np = this->the_points.size();
    uint32_t nl = this->the_lines.size();
    FILE *file = fopen(filename.c_str(),"w+");
    fprintf(file,"# vtk DataFile Version 4.1\n");
    fprintf(file,"vtk output\n");
    fprintf(file,"ASCII\n");
    fprintf(file,"DATASET POLYDATA\n");
    fprintf(file,"POINTS %u float\n",np);
    for (uint32_t i = 0; i < np; i++)
        fprintf(file,"%g %g %g\n",this->the_points[i].pos[0],this->the_points[i].pos[1],this->the_points[i].pos[2]);
    fprintf(file,"LINES %u %u\n",nl,nl*3);
    for (uint32_t i = 0; i < nl; i++)
        fprintf(file,"2 %u %u\n",this->the_lines[i].src,this->the_lines[i].dest);
    if (this->has_point_data)
    {
        fprintf(file,"POINT_DATA %u\n",np);
        fprintf(file,"SCALARS sigma float\n");
        fprintf(file,"LOOKUP_TABLE default\n");
        for (uint32_t i = 0; i < np; i++)
            fprintf(file,"%g\n",this->the_point_scalars[i]);
    }
    fclose(file);
}
