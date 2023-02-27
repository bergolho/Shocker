//
// Created by bergolho on 12/02/19.
//

#ifndef POINT_H
#define POINT_H

#include <cstdio>
#include <cstdlib>
#include <cstdbool>
#include <cstdint>
#include <cstring>
#include <vector>

class Point_Segment
{
public:
    uint32_t id;
    double pos[3];
    bool is_active;
public:
    Point_Segment () { };
    Point_Segment (const uint32_t id) { this->id = id; }
    Point_Segment (const uint32_t id, const double pos[]) { this->id = id; memcpy(this->pos,pos,sizeof(double)*3); this->is_active = false; }
    Point_Segment (Point_Segment *input) { this->id = input->id; memcpy(this->pos,input->pos,sizeof(double)*3); this->is_active = input->is_active; }
    inline void setId (const uint32_t id) { this->id = id; }
    inline void setPosition (const double pos[]) { memcpy(this->pos,pos,sizeof(double)*3); }
    inline void setActive (const bool is_active) { this->is_active = is_active; }
    Point_Segment* copy () 
    { 
        Point_Segment *result = new Point_Segment(); 
        result->id = this->id;
        memcpy(result->pos,this->pos,sizeof(double)*3);
        result->is_active = this->is_active; 
        return result; 
    }
    void print () { printf("id=%u, pos=(%g %g %g), is_active=%d\n",this->id,this->pos[0],this->pos[1],this->pos[2],static_cast<int>(this->is_active)); }
};

Point_Segment* search_point (std::vector<Point_Segment*> p_list, const uint32_t index);
void eliminate_point_from_list (std::vector<Point_Segment*> &p_list, Point_Segment *p);
void order_point_list (std::vector<Point_Segment*> &p_list);

#endif