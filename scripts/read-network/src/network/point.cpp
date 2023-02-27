#include "point.h"

void eliminate_point_from_list (std::vector<Point_Segment*> &p_list, Point_Segment *p)
{
    p_list.erase(p_list.begin() + p->id);
    order_point_list(p_list);
    delete p;
    p = NULL;
}

void order_point_list (std::vector<Point_Segment*> &p_list)
{
    uint32_t counter = 0;
    for (uint32_t i = 0; i < p_list.size(); i++)
    {
        p_list[i]->id = counter;
        counter++;
    }
}

Point_Segment* search_point (std::vector<Point_Segment*> p_list, const uint32_t index)
{
    return p_list[index];
}