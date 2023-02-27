#ifndef CLOUD_CONFIG_H
#define CLOUD_CONFIG_H

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <string>

class CloudConfig
{
public:
    std::string cloud_filename;
    std::string surface_filename;
    uint32_t rand_offset;
    double phi;
    double radius;
public:
    CloudConfig ();
    ~CloudConfig ();
    void print ();
};

#endif