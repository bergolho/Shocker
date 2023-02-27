#include "logger.h"

Logger::Logger (User_Options *options)
{
    this->iter = 0;
    write_configuration_file(options);
}

Logger::~Logger ()
{

}

void Logger::write_full_network_to_vtk (const std::vector<Node*> &n_list, const std::vector<Segment*> &s_list,\
                                    const std::string output_dir, const double his_offset, const uint32_t num_terminals)
{
    uint32_t num_nodes = n_list.size();
    uint32_t num_segments = s_list.size();
    
    std::string filename;
    std::ostringstream os;
    //os << this->params->output_dir << "/refinement_iter_" << refinement_iter << "/tree_nterm_" << num_terminals << ".vtk";
    //os << output_dir << "/tree_nterm_" << num_terminals << ".vtk";
    os << output_dir << "/tree_iter_" << this->iter << ".vtk";
    filename = os.str();

    FILE *file = fopen(filename.c_str(),"w+");
    fprintf(file,"# vtk DataFile Version 3.0\n");
    fprintf(file,"Tree\n");
    fprintf(file,"ASCII\n");
    fprintf(file,"DATASET POLYDATA\n");
    fprintf(file,"POINTS %u float\n",num_nodes);
    for (uint32_t i = 0; i < num_nodes; i++)
        fprintf(file,"%g %g %g\n",n_list[i]->pos[0],n_list[i]->pos[1],n_list[i]->pos[2]);
    fprintf(file,"LINES %u %u\n",num_segments,num_segments*3);
    for (uint32_t i = 0; i < num_segments; i++)
        fprintf(file,"2 %u %u\n",s_list[i]->src->id,s_list[i]->dest->id);
    fprintf(file,"CELL_DATA %u\n",num_segments);
    fprintf(file,"FIELD FieldData 4\n");
    fprintf(file,"LAT 1 %u float\n",num_segments);
    for (uint32_t i = 0; i < num_segments; i++)
        fprintf(file,"%g\n",s_list[i]->lat + his_offset);
    fprintf(file,"METADATA\n");
    fprintf(file,"INFORMATION 0\n\n");
    fprintf(file,"diameter 1 %u float\n",num_segments);
    for (uint32_t i = 0; i < num_segments; i++)
        fprintf(file,"%g\n",s_list[i]->diameter);
    fprintf(file,"METADATA\n");
    fprintf(file,"INFORMATION 0\n\n");
    fprintf(file,"length 1 %u float\n",num_segments);
    for (uint32_t i = 0; i < num_segments; i++)
        fprintf(file,"%g\n",s_list[i]->length);
    fprintf(file,"METADATA\n");
    fprintf(file,"INFORMATION 0\n\n");
    fprintf(file,"minTree 1 %u float\n",num_segments);
    for (uint32_t i = 0; i < num_segments; i++)
        fprintf(file,"%d\n",static_cast<int>(s_list[i]->mintree));
    fprintf(file,"METADATA\n");
    fprintf(file,"INFORMATION 0\n\n");
    fclose(file);

    this->iter++;
}

void Logger::write_full_network_monoalg_to_vtk (const std::vector<Node*> &n_list, const std::vector<Segment*> &s_list,\
                                    const std::string output_dir, const double his_offset, const uint32_t num_terminals)
{
    uint32_t num_nodes = n_list.size();
    uint32_t num_segments = s_list.size();
    
    std::string filename;
    std::ostringstream os;
    //os << this->params->output_dir << "/refinement_iter_" << refinement_iter << "/tree_nterm_" << num_terminals << ".vtk";
    os << output_dir << "/tree_nterm_monoalg" << num_terminals << ".vtk";
    filename = os.str();

    double sigmas[num_nodes];
    double root_sigma = s_list[0]->calc_conductivity();
    for (uint32_t i = 0; i < num_nodes; i++)
        sigmas[i] = root_sigma;
    
    for (uint32_t i = 0; i < s_list.size(); i++)
    {
        Segment *s = s_list[i];
        Node *dest = s->dest;
        double sigma = s->calc_conductivity();
        sigmas[dest->id] = sigma;
    }

    FILE *file = fopen(filename.c_str(),"w+");
    fprintf(file,"# vtk DataFile Version 3.0\n");
    fprintf(file,"Tree\n");
    fprintf(file,"ASCII\n");
    fprintf(file,"DATASET POLYDATA\n");
    fprintf(file,"POINTS %u float\n",num_nodes);
    for (uint32_t i = 0; i < num_nodes; i++)
        fprintf(file,"%g %g %g\n",n_list[i]->pos[0],n_list[i]->pos[1],n_list[i]->pos[2]);
    fprintf(file,"LINES %u %u\n",num_segments,num_segments*3);
    for (uint32_t i = 0; i < num_segments; i++)
        fprintf(file,"2 %u %u\n",s_list[i]->src->id,s_list[i]->dest->id);
    fprintf(file,"POINT_DATA %u\n",num_nodes);
    fprintf(file,"SCALARS sigma float\n");
    fprintf(file,"LOOKUP_TABLE default\n");
    for (uint32_t i = 0; i < num_nodes; i++)
        fprintf(file,"%g\n",sigmas[i]);
    fprintf(file,"METADATA\n");
    fprintf(file,"INFORMATION 0\n\n");
    fclose(file);
}

void Logger::write_minimum_network_to_vtk (const std::vector<Node*> &n_list, const std::vector<Segment*> &s_list,\
                                    const std::string output_dir, const double his_offset, const uint32_t num_terminals)
{
    uint32_t num_nodes = n_list.size();
    uint32_t num_segments = s_list.size();
    
    std::string filename;
    std::ostringstream os;
    os << output_dir << "/minimum_network.vtk";
    filename = os.str();

    FILE *file = fopen(filename.c_str(),"w+");
    fprintf(file,"# vtk DataFile Version 3.0\n");
    fprintf(file,"Tree\n");
    fprintf(file,"ASCII\n");
    fprintf(file,"DATASET POLYDATA\n");
    fprintf(file,"POINTS %u float\n",num_nodes);
    for (uint32_t i = 0; i < num_nodes; i++)
        fprintf(file,"%g %g %g\n",n_list[i]->pos[0],n_list[i]->pos[1],n_list[i]->pos[2]);
    fprintf(file,"LINES %u %u\n",num_segments,num_segments*3);
    for (uint32_t i = 0; i < num_segments; i++)
        fprintf(file,"2 %u %u\n",s_list[i]->src->id,s_list[i]->dest->id);
    fprintf(file,"CELL_DATA %u\n",num_segments);
    fprintf(file,"FIELD FieldData 4\n");
    fprintf(file,"LAT 1 %u float\n",num_segments);
    for (uint32_t i = 0; i < num_segments; i++)
        fprintf(file,"%g\n",s_list[i]->lat + his_offset);
    fprintf(file,"METADATA\n");
    fprintf(file,"INFORMATION 0\n\n");
    fprintf(file,"diameter 1 %u float\n",num_segments);
    for (uint32_t i = 0; i < num_segments; i++)
        fprintf(file,"%g\n",s_list[i]->diameter);
    fprintf(file,"METADATA\n");
    fprintf(file,"INFORMATION 0\n\n");
    fprintf(file,"length 1 %u float\n",num_segments);
    for (uint32_t i = 0; i < num_segments; i++)
        fprintf(file,"%g\n",s_list[i]->length);
    fprintf(file,"METADATA\n");
    fprintf(file,"INFORMATION 0\n\n");
    fprintf(file,"minTree 1 %u float\n",num_segments);
    for (uint32_t i = 0; i < num_segments; i++)
        fprintf(file,"%d\n",static_cast<int>(s_list[i]->mintree));
    fprintf(file,"METADATA\n");
    fprintf(file,"INFORMATION 0\n\n");
    fclose(file);
}

void Logger::write_minimum_network_monoalg_to_vtk (const std::vector<Node*> &n_list, const std::vector<Segment*> &s_list,\
                                    const std::string output_dir, const double his_offset, const uint32_t num_terminals)
{
    uint32_t num_nodes = n_list.size();
    uint32_t num_segments = s_list.size();
    
    std::string filename;
    std::ostringstream os;
    os << output_dir << "/minimum_network_monoalg.vtk";
    filename = os.str();

    double sigmas[num_nodes];
    double root_sigma = s_list[0]->calc_conductivity();
    for (uint32_t i = 0; i < num_nodes; i++)
        sigmas[i] = root_sigma;

    for (uint32_t i = 0; i < num_segments; i++)
    {
        Segment *s = s_list[i];
        Node *dest = s->dest;
        double sigma = s->calc_conductivity();
        sigmas[dest->id] = sigma;
    }

    FILE *file = fopen(filename.c_str(),"w+");
    fprintf(file,"# vtk DataFile Version 3.0\n");
    fprintf(file,"Tree\n");
    fprintf(file,"ASCII\n");
    fprintf(file,"DATASET POLYDATA\n");
    fprintf(file,"POINTS %u float\n",num_nodes);
    for (uint32_t i = 0; i < num_nodes; i++)
        fprintf(file,"%g %g %g\n",n_list[i]->pos[0],n_list[i]->pos[1],n_list[i]->pos[2]);
    fprintf(file,"LINES %u %u\n",num_segments,num_segments*3);
    for (uint32_t i = 0; i < num_segments; i++)
        fprintf(file,"2 %u %u\n",s_list[i]->src->id,s_list[i]->dest->id);
    fprintf(file,"POINT_DATA %u\n",num_nodes);
    fprintf(file,"SCALARS sigma float\n");
    fprintf(file,"LOOKUP_TABLE default\n");
    for (uint32_t i = 0; i < num_nodes; i++)
        fprintf(file,"%g\n",sigmas[i]);
    fprintf(file,"METADATA\n");
    fprintf(file,"INFORMATION 0\n\n");
    fclose(file);
}

void Logger::write_simulation_time (const std::string output_dir, const long res_time, const uint32_t inactive_points_time, const uint32_t active_points_time,\
                                const uint32_t inactive_search_time, const uint32_t active_search_time,\
                                const uint32_t inactive_eval_time, const uint32_t active_eval_time,\
                                const uint32_t main_loop_time, const uint32_t force_connection_no_lat_error_tolerance_time, const uint32_t force_connection_no_distance_criterion_time,\
                                const uint32_t force_connection_no_distance_criterion_geodesic_time, const uint32_t force_connection_no_distance_criterion_line_time,\
                                const uint32_t geodesic_pathway_time,\
                                const double conv_rate)
{
    std::string filename;
    std::ostringstream os;
    os << output_dir << "/simulation_time.txt";
    filename = os.str();

    FILE *file = fopen(filename.c_str(),"w+");
    fprintf(file,"[INFO] Resolution time = %lu μs (%lf min)\n",res_time,res_time/conv_rate);
    fprintf(file,"[INFO] Total inactive points time = %u μs (%lf min)\n",inactive_points_time,inactive_points_time/conv_rate);
    fprintf(file,"[INFO] Total active points time = %u μs (%lf min)\n",active_points_time,active_points_time/conv_rate);
    fprintf(file,"[INFO] Total inactive points search time = %u μs (%lf min)\n",inactive_search_time,inactive_search_time/conv_rate);
    fprintf(file,"[INFO] Total active points search time = %u μs (%lf min)\n",active_search_time,active_search_time/conv_rate);
    fprintf(file,"[INFO] Total inactive points evaluation time = %u μs (%lf min)\n",inactive_eval_time,inactive_eval_time/conv_rate);
    fprintf(file,"[INFO] Total active points evaluation time = %u μs (%lf min)\n",active_eval_time,active_eval_time/conv_rate);
    fprintf(file,"[INFO] Total geodesic pathway time = %u μs (%lf min)\n",geodesic_pathway_time,geodesic_pathway_time/conv_rate);
    fprintf(file,"[INFO] Total main loop time = %u μs (%lf min)\n",main_loop_time,main_loop_time/conv_rate);
    fprintf(file,"[INFO] Total force connection without LAT error tolerance time = %u μs (%lf min)\n",force_connection_no_lat_error_tolerance_time,force_connection_no_lat_error_tolerance_time/conv_rate);
    fprintf(file,"[INFO] Total force connection without distance criterion time = %u μs (%lf min)\n",force_connection_no_distance_criterion_time,force_connection_no_distance_criterion_time/conv_rate);
    fprintf(file,"[INFO] Total force connection without distance criterion time (geodesic) = %u μs (%lf min)\n",force_connection_no_distance_criterion_geodesic_time,force_connection_no_distance_criterion_geodesic_time/conv_rate);
    fprintf(file,"[INFO] Total force connection without distance criterion time (line) = %u μs (%lf min)\n",force_connection_no_distance_criterion_line_time,force_connection_no_distance_criterion_line_time/conv_rate);
    
    fclose(file);
}

void Logger::write_configuration_file (User_Options *options)
{
    std::string filename;
    std::ostringstream os;
    os << options->output_dir << "/configuration_file.ini";
    filename = os.str();

    printf("[logger] Configuration file to reproduce this simulation:> %s\n",filename.c_str());
    FILE *file = fopen(filename.c_str(),"w+");
    fprintf(file,"[main]\n");
    fprintf(file,"root_pos = [%lf, %lf, %lf]\n",options->root_pos[0],options->root_pos[1],options->root_pos[2]);
    fprintf(file,"seed = %u\n",options->seed);
    fprintf(file,"start_diameter = %lf\n",options->start_diameter);
    fprintf(file,"l_d = %lf\n",options->characteristic_length);
    fprintf(file,"his_offset = %lf\n",options->his_offset);
    fprintf(file,"N_p = %u\n",options->np);
    fprintf(file,"N_a = %u\n",options->na);
    if (options->use_initial_network)
    {
        fprintf(file,"use_initial_network = true\n");
        fprintf(file,"initial_network_filename = %s\n",options->initial_network_filename.c_str());
    }
    fprintf(file,"\n");
    fprintf(file,"[save_network]\n");
    fprintf(file,"output_dir = %s\n",options->output_dir.c_str());
    fprintf(file,"\n");
    fprintf(file,"[cloud_points]\n");
    fprintf(file,"phi = %lf\n",options->cloud_config->phi);
    fprintf(file,"rand_offset = %u\n",options->cloud_config->rand_offset);
    fprintf(file,"use_cloud_point = true\n");
    fprintf(file,"cloud_points_filename = %s\n",options->cloud_config->cloud_filename.c_str());
    fprintf(file,"surface_filename = %s\n",options->cloud_config->surface_filename.c_str());
    fprintf(file,"\n");
    fprintf(file,"[pmj]\n");
    if (options->use_pmj_location)
    {
        fprintf(file,"use_pmj_location = true\n");
        fprintf(file,"pmj_location_filename = %s\n",options->pmj_config->location_filename.c_str());
        fprintf(file,"pmj_connection_rate = %u\n",options->pmj_config->connection_rate);
        fprintf(file,"lat_error_tolerance = %lf\n",options->pmj_config->lat_error_tolerance);
    }
    else
    {
        fprintf(file,"use_pmj_location = false\n");
    }
    fprintf(file,"\n");
    fprintf(file,"[cost_function]\n");
    fprintf(file,"function_name = custom_function\n");
    fprintf(file,"beta = 1.0\n");
    fprintf(file,"alpha = 0.0\n");
    fprintf(file,"min_degrees_limit = 10\n");
    fprintf(file,"max_degrees_limit = 90\n");
    fprintf(file,"min_segment_limit = 100\n");
    fprintf(file,"max_segment_limit = 10000\n");
    fclose(file);
}

void Logger::write_lmin_to_file (const std::string output_dir, const uint32_t num_terminal, const double l_min)
{
    std::string filename;
    std::ostringstream os;
    
    os << output_dir << "/l_min.txt";
    filename = os.str();

    FILE *file = fopen(filename.c_str(),"a");
    fprintf(file,"%u %g\n",num_terminal,l_min);
    fclose(file);
}