//
// Created by bergolho on 12/02/19.
//

#ifndef SEGMENT_H
#define SEGMENT_H

#include <cstdio>
#include <cstdlib>
#include <cstdbool>
#include <cstdint>
#include <cmath>

#include "constants.h"
#include "point.h"

class Segment
{
public:
    uint32_t id;
    double diameter;
    double length;
    double lat;

    bool mintree;

    Point_Segment *src;
    Point_Segment *dest;

    Segment *left;
    Segment *right;
    Segment *parent;
public:
    Segment () { }
    Segment (const uint32_t id, Point_Segment *src, Point_Segment *dest, Segment *left, Segment *right, Segment *parent)
    {
        this->id = id;
        this->src = src;
        this->length = sqrt( powf(src->pos[0]-dest->pos[0],2) + powf(src->pos[1]-dest->pos[1],2) + powf(src->pos[2]-dest->pos[2],2) );
        this->dest = dest;
        this->left = left;
        this->right = right;
        this->parent = parent;
        this->diameter = 0.0;
        this->lat = 0.0;
        this->mintree = false;
    }
    ~Segment ()
    {
        this->parent = NULL;
        this->left = NULL;
        this->right = NULL;
        this->src = NULL;
        this->dest = NULL;
    }
    void calc_middle_point (double pos[])
    {
        Point_Segment *src = this->src;
        Point_Segment *dest = this->dest;

        pos[0] = (src->pos[0] + dest->pos[0]) / 2.0;
        pos[1] = (src->pos[1] + dest->pos[1]) / 2.0;
        pos[2] = (src->pos[2] + dest->pos[2]) / 2.0;        
    }
    void calc_unitary_vector (double u[])
    {
        if (this->length < 1.0e-08) this->length = 1.0e-08;

        u[0] = (this->dest->pos[0] - this->src->pos[0]) / this->length;
        u[1] = (this->dest->pos[1] - this->src->pos[1]) / this->length;
        u[2] = (this->dest->pos[2] - this->src->pos[2]) / this->length;
    }
    double calc_diameter ()
    {
        return this->diameter;
    }
    double calc_length ()
    {
        return sqrt( pow(src->pos[0]-dest->pos[0],2) + pow(src->pos[1]-dest->pos[1],2) + pow(src->pos[2]-dest->pos[2],2) );
    }
    double calc_propagation_velocity ()
    {
        const double G = 7.9;
        const double Cf = 3.4;
        const double tauf = 0.1;

        double d = this->diameter;
    
        // Output in {m/s}
        return pow( (G*d)/(4.0*Cf*tauf) , 0.5 ) * 0.1;
    }
    double calc_local_activation_time ()
    {
        double length = this->length;
        double cv = calc_propagation_velocity()*M_S_TO_UM_MS;
        return length*M_TO_UM / cv;
    }
    double calc_pathway_length ()
    {
        double result = 0.0;
        Segment *tmp = this;
        while (tmp != NULL)
        {
            result += tmp->length;
            tmp = tmp->parent;
        }
        return result;
    }
    double calc_terminal_local_activation_time ()
    {
        double result = 0.0;
        Segment *tmp = this;
        while (tmp != NULL)
        {
            double cv = tmp->calc_propagation_velocity()*M_S_TO_UM_MS;    // {m/s}->{um/ms}
            double dist = tmp->length;                            
            double lat = dist/cv;

            result += lat;
            
            tmp = tmp->parent;
        }

        // Return the LAT in {ms}
        return result; 
    }
    uint32_t calc_level ()
    {
        uint32_t result = 0;
        Segment *tmp = this;
        while (tmp != NULL)
        {
            result++;
            tmp = tmp->parent;
        }
        return result;
    }
    bool is_terminal ()
    {
        return (this->left == NULL && this->right == NULL) ? true : false;
    }
    bool is_bifurcation ()
    {
        return (this->left != NULL && this->right != NULL) ? true : false;
    }
    bool is_inside_region (const double center[], const double radius)
    {
        double dist = sqrt( pow(dest->pos[0]-center[0],2) + pow(dest->pos[1]-center[1],2) + pow(dest->pos[2]-center[2],2) ); 
        return (dist < radius) ? true : false;
    }
    bool is_mintree ()
    {
        return this->mintree;
    }
    void print ()
    {
        printf("Segment %u (%u,%u) (%g %g %g) - (%g %g %g) -- DIAMETER = %g um -- LAT = %g ms\n\n",this->id,src->id,dest->id,\
                                                                            src->pos[0],src->pos[1],src->pos[2],dest->pos[0],dest->pos[1],dest->pos[2],\
                                                                            this->diameter,this->lat);
        if (this->parent == NULL)
            printf("\tPARENT = NIL\n");
        else
            printf("\tPARENT = %u\n",this->parent->id);
        if (this->left == NULL)
            printf("\tLEFT = NIL\n");
        else
            printf("\tLEFT = %u\n",this->left->id);
        if (this->right == NULL)
            printf("\tRIGHT = NIL\n");
        else
            printf("\tRIGHT = %u\n",this->right->id);
    }
};

void eliminate_segment_from_list (std::vector<Segment*> &s_list, Segment *s);
void order_segment_list (std::vector<Segment*> &s_list);

#endif