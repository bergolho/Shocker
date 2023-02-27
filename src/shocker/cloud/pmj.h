#ifndef PMJ_DATA_H
#define PMJ_DATA_H

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>

#include <vector>
#include <string>
#include <algorithm>

#include "../../utils/utils.h"
#include "../../options/pmj_config.h"

class PMJPoint
{
public:
    uint32_t id;                                            // Index of a PMJ
    bool connected;                                         // Flag which tells if a PMJ is connected or not
    double pos[3];                                          // Coordinates (x,y,z) of the PMJ
    double ref_value;                                       // Reference value to be matched
    double aprox_value;                                     // Aproximated value
    double error;                                           // Error value
public:
    PMJPoint ();
    PMJPoint (const uint32_t id, const double pos[],\
        const double ref_value = 0.0, const double aprox_value = 0.0, const double error = __DBL_MAX__,\
        const bool connected = false);
    ~PMJPoint ();
    PMJPoint* copy ();
    void print ();
};

bool compareByLAT (PMJPoint *a, PMJPoint *b);

#endif