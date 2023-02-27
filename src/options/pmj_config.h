#ifndef PMJ_CONFIG_H
#define PMJ_CONFIG_H

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <string>

class PMJConfig
{
public:
    bool using_pmj;
    uint32_t connection_rate;
    double lat_error_tolerance;
    std::string location_filename;
public:
    PMJConfig ();
    ~PMJConfig ();
    void print ();
};

#endif