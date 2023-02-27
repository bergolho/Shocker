#include "segment.h"

void eliminate_segment_from_list (std::vector<Segment*> &s_list, Segment *s)
{
    s_list.erase(s_list.begin() + s->id);
    order_segment_list(s_list);
    delete s;
    s = NULL;
}

void order_segment_list (std::vector<Segment*> &s_list)
{
    uint32_t counter = 0;
    for (uint32_t i = 0; i < s_list.size(); i++)
    {
        s_list[i]->id = counter;
        counter++;
    }
}