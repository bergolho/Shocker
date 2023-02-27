#include "cloud_config.h"

CloudConfig::CloudConfig ()
{
    this->cloud_filename = "";
    this->surface_filename = "";
    this->rand_offset = 2;
    this->phi = 0.05;
    this->radius = 5000.0;
}

CloudConfig::~CloudConfig ()
{

}

void CloudConfig::print ()
{
    printf("[cloud_config] cloud_points_filename = %s\n",this->cloud_filename.c_str());
    printf("[cloud_config] surface_filename = %s\n",this->surface_filename.c_str());
    printf("[cloud_config] phi = %g\n",this->phi);
    printf("[cloud_config] radius = %g\n",this->radius);
    printf("[cloud_config] rand_offset = %u\n",this->rand_offset);
}