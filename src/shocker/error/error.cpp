#include "error.h"

Error::Error ()
{
    this->counter_main_loop_connections = 0;
    this->counter_no_lat_error_tolerance_connections = 0;
    this->counter_no_distance_criterion_connections = 0;
    this->counter_no_distance_criterion_connections_geodesic = 0;
    this->counter_no_distance_criterion_connections_line = 0;

    this->max_lat_error = __DBL_MIN__;
    this->min_max_aprox_lat[0] = __DBL_MAX__; this->min_max_aprox_lat[1] = __DBL_MIN__;
    this->min_max_ref_lat[0] = __DBL_MAX__; this->min_max_ref_lat[1] = __DBL_MIN__;
    this->min_max_term_lat[0] = __DBL_MAX__; this->min_max_term_lat[1] = __DBL_MIN__;

    this->rmse = __DBL_MAX__;
    this->rrmse = __DBL_MAX__;
    this->epsilon_2ms = __UINT32_MAX__;
    this->epsilon_5ms = __UINT32_MAX__;
}

Error::~Error ()
{

}

void Error::calculate_electric_error (PMJ_Data *pmj_data)
{
    // Compute min/max LAT values from the active PMJ's that are connected to the network
    // Compute max LAT error from the active PMJ's that are connected to the network
    uint32_t n = pmj_data->pmjs.size();
    for (uint32_t i = 0; i < n; i++)
    {
        if (pmj_data->pmjs[i]->connected)
        {
            double aprox_lat = pmj_data->pmjs[i]->aprox_value;
            if (aprox_lat < this->min_max_aprox_lat[0]) this->min_max_aprox_lat[0] = aprox_lat;
            if (aprox_lat > this->min_max_aprox_lat[1]) this->min_max_aprox_lat[1] = aprox_lat;

            double ref_lat = pmj_data->pmjs[i]->ref_value;
            if (ref_lat < this->min_max_ref_lat[0]) this->min_max_ref_lat[0] = ref_lat;
            if (ref_lat > this->min_max_ref_lat[1]) this->min_max_ref_lat[1] = ref_lat;

            double error = fabs(pmj_data->pmjs[i]->error);
            if (error > this->max_lat_error) this->max_lat_error = error;
        }   
    }

    // Compute the RMSE and RRMSE from the active PMJ's
    uint32_t counter = 0;
    double sum_num = 0.0;
    double sum_den = 0.0;
    for (uint32_t i = 0; i < n; i++)
    {
        if (pmj_data->pmjs[i]->connected)
        {
            double ref_value = pmj_data->pmjs[i]->ref_value;
            double error = fabs(pmj_data->pmjs[i]->error);

            sum_num += powf(error,2);
            sum_den += powf(ref_value,2);
            counter++;
        }
    }    
    double l2_norm = sqrt(sum_den);
    this->rmse = sqrt(sum_num/(double)counter);
    this->rrmse = sqrt(sum_num/sum_den);

    // Compute the number of active PMJ's that have an error less than a certain threashold
    uint32_t counter_less_2ms = 0;
    uint32_t counter_less_5ms = 0;
    for (uint32_t i = 0; i < n; i++)
    {
        if (pmj_data->pmjs[i]->connected)
        {
            double error = fabs(pmj_data->pmjs[i]->error);

            if (error < 2.0)    counter_less_2ms++;
            if (error < 5.0)    counter_less_5ms++;
        }
        
    }
    this->epsilon_2ms = (double)counter_less_2ms / (double)counter;
    this->epsilon_5ms = (double)counter_less_5ms / (double)counter;
}

void Error::update_min_max_terminal_lat (std::vector<Segment*> s_list)
{
    for (uint32_t i = 0; i < s_list.size(); i++)
    {
        Segment *cur_segment = s_list[i];
        if (cur_segment->is_terminal())
        {
            double lat = cur_segment->calc_terminal_local_activation_time();
            if (lat < this->min_max_term_lat[0]) this->min_max_term_lat[0] = lat;
            if (lat > this->min_max_term_lat[1]) this->min_max_term_lat[1] = lat;
        }
    }
    std::cout << get_color_code(GRAY) << "[error] Terminals --> min.LAT = " << this->min_max_term_lat[0] << " || max.LAT = " << this->min_max_term_lat[1] << std::endl;
}

Error* Error::copy ()
{
    Error *result = new Error();

    result->max_lat_error = this->max_lat_error;
    result->counter_main_loop_connections = this->counter_main_loop_connections;
    result->counter_no_lat_error_tolerance_connections = this->counter_no_lat_error_tolerance_connections;
    result->counter_no_distance_criterion_connections = this->counter_no_distance_criterion_connections;
    result->counter_no_distance_criterion_connections_geodesic = this->counter_no_distance_criterion_connections_geodesic;
    result->counter_no_distance_criterion_connections_line = this->counter_no_distance_criterion_connections_line;
    memcpy(result->min_max_aprox_lat,this->min_max_aprox_lat,sizeof(double)*2);
    memcpy(result->min_max_ref_lat,this->min_max_ref_lat,sizeof(double)*2);
    memcpy(result->min_max_term_lat,this->min_max_term_lat,sizeof(double)*2);
    
    return result;
}

void Error::concatenate (Error *input)
{
    this->counter_main_loop_connections += input->counter_main_loop_connections;
    this->counter_no_lat_error_tolerance_connections += input->counter_no_lat_error_tolerance_connections;
    this->counter_no_distance_criterion_connections += input->counter_no_distance_criterion_connections;
    this->counter_no_distance_criterion_connections_geodesic += input->counter_no_distance_criterion_connections_geodesic;
    this->counter_no_distance_criterion_connections_line += input->counter_no_distance_criterion_connections_line;

    this->min_max_aprox_lat[0] = (input->min_max_aprox_lat[0] < this->min_max_aprox_lat[0]) ? input->min_max_aprox_lat[0] : this->min_max_aprox_lat[0];
    this->min_max_aprox_lat[1] = (input->min_max_aprox_lat[1] > this->min_max_aprox_lat[1]) ? input->min_max_aprox_lat[1] : this->min_max_aprox_lat[1];
    this->min_max_ref_lat[0] = (input->min_max_ref_lat[0] < this->min_max_ref_lat[0]) ? input->min_max_ref_lat[0] : this->min_max_ref_lat[0];
    this->min_max_ref_lat[1] = (input->min_max_ref_lat[1] > this->min_max_ref_lat[1]) ? input->min_max_ref_lat[1] : this->min_max_ref_lat[1];
    this->max_lat_error = (input->max_lat_error > this->max_lat_error) ? input->max_lat_error : this->max_lat_error;
}

void Error::reset ()
{
    this->counter_main_loop_connections = 0;
    this->counter_no_lat_error_tolerance_connections = 0;
    this->counter_no_distance_criterion_connections = 0;
    this->counter_no_distance_criterion_connections_geodesic = 0;
    this->counter_no_distance_criterion_connections_line = 0;

    this->max_lat_error = __DBL_MIN__;
    this->min_max_aprox_lat[0] = __DBL_MAX__; this->min_max_aprox_lat[1] = __DBL_MIN__;
    this->min_max_ref_lat[0] = __DBL_MAX__; this->min_max_ref_lat[1] = __DBL_MIN__;
    this->min_max_term_lat[0] = __DBL_MAX__; this->min_max_term_lat[1] = __DBL_MIN__;

    this->rmse = __DBL_MAX__;
    this->rrmse = __DBL_MAX__;
    this->epsilon_2ms = __UINT32_MAX__;
    this->epsilon_5ms = __UINT32_MAX__;
}

void Error::print ()
{
    printf("[error] Min/Max terminal LAT = [ %g %g ]\n",this->min_max_term_lat[0],this->min_max_term_lat[1]);
    printf("[error] Min/Max aproximation LAT = [ %g %g ]\n",this->min_max_aprox_lat[0],this->min_max_aprox_lat[1]);
    printf("[error] Min/Max reference LAT = [ %g %g ]\n",this->min_max_ref_lat[0],this->min_max_ref_lat[1]);
    printf("[error] Max error LAT = %g\n",this->max_lat_error);
    printf("[error] RMSE = %g\n",this->rmse);
    printf("[error] RRMSE = %g\n",this->rrmse);
    printf("[error] Epsilon < 2ms = %g\n",this->epsilon_2ms);
    printf("[error] Epsilon < 5ms = %g\n",this->epsilon_5ms);
}