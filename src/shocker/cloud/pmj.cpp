#include "pmj.h"

PMJPoint::PMJPoint ()
{
    this->id = 0;
    this->connected = false;
    this->ref_value = 0.0;
    this->aprox_value = 0.0;
    this->error = __DBL_MAX__;
}

PMJPoint::PMJPoint (const uint32_t id, const double pos[], const double ref_value, const double aprox_value, const double error, const bool connected)
{
    this->id = id;
    this->connected = connected;
    this->ref_value = ref_value;
    this->aprox_value = aprox_value;
    this->error = error;
    memcpy(this->pos,pos,sizeof(double)*3);
}

PMJPoint::~PMJPoint ()
{

}

PMJPoint* PMJPoint::copy ()
{
    PMJPoint *result = new PMJPoint(this->id,this->pos,this->ref_value,this->aprox_value,this->error,this->connected);
    return result;
}

void PMJPoint::print ()
{
    printf("[pmj] id = %u || pos = (%g %g %g) || Ref = %g || Aprox = %g || Error = %g || Connected = %d\n",\
            this->id,this->pos[0],this->pos[1],this->pos[2],this->ref_value,this->aprox_value,this->error,(int)this->connected);
}

bool compareByLAT (PMJPoint *a, PMJPoint *b)
{
    return (a->ref_value < b->ref_value);
}