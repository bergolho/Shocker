#include "pmj_config.h"

PMJConfig::PMJConfig ()
{
    this->using_pmj = false;
    this->lat_error_tolerance = 2.0;
    this->connection_rate = __UINT32_MAX__;
    this->location_filename = "";
}

PMJConfig::~PMJConfig ()
{

}

void PMJConfig::print ()
{
    printf("LAT tolerance error = %g ms\n",this->lat_error_tolerance);
    printf("PMJ connection rate = %u\n",this->connection_rate);
    printf("PMJ location filename = \"%s\"\n",this->location_filename.c_str());
}