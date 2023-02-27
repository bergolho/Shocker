//
// Created by bergolho on 18/08/21.
//

#ifndef CLOUD_POINT_H
#define CLOUD_POINT_H

#include <cstdio>
#include <cstdlib>
#include <cstdbool>
#include <cstdint>
#include <cstring>
#include <cmath>

#include <vector>
#include <set>
#include <string>
#include <algorithm>

#include "pmj.h"

class CloudPoint
{
public:
    bool taken;
    uint32_t id;
    double pos[3];
    double lat;
    PMJPoint *closest_pmj;
public:
    CloudPoint ();
    CloudPoint (const uint32_t id, const double pos[], const double lat, const bool taken);
    CloudPoint (CloudPoint *in);
    ~CloudPoint ();
    inline void setTaken (const bool taken) { this->taken = taken; }
    inline void setId (const uint32_t id) { this->id = id; }
    inline void setPosition (const double pos[]) { memcpy(this->pos,pos,sizeof(double)*3); }
    inline void setLAT (const double lat) { this->lat = lat; }
    void print ();
};

#endif