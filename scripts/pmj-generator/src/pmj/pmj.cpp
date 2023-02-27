#include "pmj.h"

PMJ::PMJ ()
{

}

PMJ::PMJ (const uint32_t id, const uint32_t term_id, const double pos[],\
         const double pk_lat, const double tiss_lat, const double delay,\
         const double pk_dist, const double pk_cv)
{
    this->id = id;
    this->term_id = term_id;
    memcpy(this->pos,pos,sizeof(double)*3);
    this->pk_lat = pk_lat;
    this->tiss_lat = tiss_lat;
    this->delay = delay;
    this->pk_dist = pk_dist;
    this->pk_cv = pk_cv;
    this->is_active = (delay == 0) ? false : true;
}

void PMJ::print ()
{
    printf("|| term_id = %u || id = %u || pos = (%g %g %g) || pk_LAT = %g || tiss_LAT = %g || pmj_delay = %g || pk_Dist = %g || pk_CV = %g ||\n",\
            this->term_id,this->id,this->pos[0],this->pos[1],this->pos[2],this->pk_lat,this->tiss_lat,this->delay,this->pk_dist,this->pk_cv);
}

void get_closest_points_inside_sphere (std::vector<Point> endo_points, std::vector<PMJ> pmj_points, std::vector<Point> &out, const double percentage, const double his_offset)
{
    const uint32_t nc = endo_points.size();
    std::vector<bool> taken;
    taken.assign(nc,false);

    for (uint32_t i = 0; i < pmj_points.size(); i++)
    {
        print_progress_bar(i,pmj_points.size());
        
        PMJ cur_pmj = pmj_points[i];
        if (cur_pmj.is_active)
        {
            for (uint32_t j = 0; j < nc; j++)
            {
                double dist = calc_norm(cur_pmj.pos[0],cur_pmj.pos[1],cur_pmj.pos[2],endo_points[j].x,endo_points[j].y,endo_points[j].z);
                if (dist < PMJ_REGION_RADIUS && !taken[j])
                {
                    taken[j] = true;

                    // Subtract the PMJ delay and the His-Bundle time from the endocardium point LAT
                    endo_points[j].value -= (cur_pmj.delay+his_offset);
                }
            }
        }
    }

    std::vector<Point> points;
    for (uint32_t i = 0; i < nc; i++)
    {
        if (taken[i])
        {
            Point p(i,endo_points[i].x,endo_points[i].y,endo_points[i].z,endo_points[i].value);
            points.push_back(p);
        }
    }

    uint32_t np = points.size();
    uint32_t closest_np = np*percentage;
    
    taken.clear();
    taken.assign(np,false);

    uint32_t counter = 0;
    while (counter < closest_np)
    {
        uint32_t id = rand() % np;
        if (!taken[id])
        {
            out.push_back(points[id]);
            taken[id] = true;
            counter++;
        }
    }
}