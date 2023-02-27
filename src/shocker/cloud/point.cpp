#include "point.h"

CloudPoint::CloudPoint ()
{
    this->taken = false;
    this->id = 0;
    memset(this->pos,0,sizeof(double)*3);
    this->lat = 0.0;
    this->closest_pmj = nullptr;
}

CloudPoint::CloudPoint (const uint32_t id, const double pos[], const double lat, const bool taken)
{
    this->id = id;
    memcpy(this->pos,pos,sizeof(double)*3);
    this->lat = lat;
    this->taken = taken;
    this->closest_pmj = nullptr;
}

CloudPoint::CloudPoint (CloudPoint *in)
{
    this->id = in->id;
    memcpy(this->pos,in->pos,sizeof(double)*3);
    this->lat = in->lat;
    this->taken = in->taken;
    this->closest_pmj = in->closest_pmj;
}

CloudPoint::~CloudPoint ()
{

}

void CloudPoint::print ()
{
    printf("id=%u, pos=(%g %g %g), LAT=%g, Taken=%d\n",this->id,this->pos[0],this->pos[1],this->pos[2],this->lat,static_cast<int>(this->taken));
}
