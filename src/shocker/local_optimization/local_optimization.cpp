#include "local_optimization.h"

BifurcationPoint::BifurcationPoint (const uint32_t id, const double pos[])
{
    this->id = id;
    memcpy(this->pos,pos,sizeof(double)*3);
}

BifurcationPoint::~BifurcationPoint ()
{

}

void BifurcationPoint::print ()
{
    printf("id = %u, pos = (%g %g %g)\n",this->id,this->pos[0],this->pos[1],this->pos[2]);
}

LocalOptimization::LocalOptimization (LocalOptimizationConfig *local_opt_config)
{
    // Initialize the bifurcation points
    uint32_t np = (NE+1)*(2+NE)/2.0;
    double pos[3] = {0,0,0};
    for (uint32_t i = 0; i < np; i++)
    {
        BifurcationPoint b_point(i,pos);
        this->biff_points.push_back(b_point);
    }
}

LocalOptimization::~LocalOptimization ()
{

}

void LocalOptimization::update_bifurcation_points (Segment *iconn, Segment *ibiff, Segment *inew)
{
    double delta_e = 1.0 / NE;

    // Get a reference to the 3 vertices surrounding the triangle area
    Node *g1 = ibiff->src;
    Node *g2 = iconn->dest;
    Node *g3 = inew->dest;

    // Build the bifurcation points for the local optimization based on the triangle composed by
    // 'iconn', 'ibiff' and 'inew'
    uint32_t counter = 0;
    for (uint32_t i = 0; i <= NE; i++)
    {
        for (uint32_t j = 0; j <= NE-i; j++)
        {
            double epsilon = i*delta_e;
            double eta = j*delta_e;

            // Build the phi array
            double phi[3];
            phi[0] = 1.0 - epsilon - eta;
            phi[1] = epsilon;
            phi[2] = eta;

            double pos[3];
            pos[0] = (phi[0]*g1->pos[0]) + (phi[1]*g2->pos[0]) + (phi[2]*g3->pos[0]);
            pos[1] = (phi[0]*g1->pos[1]) + (phi[1]*g2->pos[1]) + (phi[2]*g3->pos[1]);
            pos[2] = (phi[0]*g1->pos[2]) + (phi[1]*g2->pos[2]) + (phi[2]*g3->pos[2]);

            this->biff_points[counter].setId(counter);
            this->biff_points[counter].setPosition(pos);
            counter++;
        }
    }
}

void LocalOptimization::print ()
{
    printf("[local_optimization] ne = %g\n",NE);
    for (uint32_t i = 0; i < this->biff_points.size(); i++)
        this->biff_points[i].print();
}

void LocalOptimization::write (std::string output_dir, const double translate[])
{
    uint32_t np = this->biff_points.size();
    std::string name = output_dir + "/bifurcation_plane_um.vtk";
    FILE *file = fopen(name.c_str(),"w+");
    fprintf(file,"# vtk DataFile Version 4.1\n");
    fprintf(file,"vtk output\n");
    fprintf(file,"ASCII\n");
    fprintf(file,"DATASET POLYDATA\n");
    fprintf(file,"POINTS %u float\n",np);
    for (uint32_t i = 0; i < np; i++)
    {
        double pos[3];
        for (uint32_t j = 0; j < 3; j++)
            pos[j] = (this->biff_points[i].pos[j] + translate[j]);
        fprintf(file,"%g %g %g\n",pos[0],pos[1],pos[2]);
    }
    fprintf(file,"VERTICES %u %u\n",np,np*2);
    for (uint32_t i = 0; i < np; i++)
        fprintf(file,"1 %u\n",i);
    fclose(file);
}