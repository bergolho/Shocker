#include "config.h"

Config::Config () { }

Config::Config (User_Options *options)
{
    set_parameters(options);
    set_save_network();
    set_cost_function(options->cost_function_config);
}

Config::~Config ()
{
    if (this->cost_fn)
        delete this->cost_fn;
    fclose(this->log_file);
}

void Config::set_parameters (User_Options *options)
{
    this->num_terminals = 0;
    this->seed = options->seed;
    
    srand(this->seed);

    memcpy(this->root_pos,options->root_pos,sizeof(double)*3);
    
    this->l_d = options->characteristic_length;
    this->his_offset = options->his_offset;
    this->start_diameter = options->start_diameter;
    this->np = options->np;
    this->na = options->na;
    
    this->using_cloud_points = options->use_cloud_points;
    this->using_pmj_location = options->use_pmj_location;
    this->using_initial_network = options->use_initial_network;

    this->output_dir = options->output_dir;
    this->initial_network_name = options->initial_network_filename;
    this->cost_function_name = options->cost_function_config->function_name;
}

void Config::set_cost_function (CostFunctionConfig *cost_fn_config)
{
    this->cost_fn = new CostFunction(cost_fn_config);
}

void Config::set_save_network ()
{
    create_directory(this->output_dir.c_str());
    printf("[config] Output directory = %s\n",this->output_dir.c_str());

    std::string log_filename = this->output_dir + "/output.log";
    this->log_file = fopen(log_filename.c_str(),"w+");
    if (!this->log_file)
    {
        fprintf(stderr,"[config] ERROR! Opening logfile!\n");
        exit(EXIT_FAILURE);
    }
    else
    {
        printf("[config] Log file will be saved at '%s'\n",log_filename.c_str());
    }    
}

void Config::set_save_network (std::string output_dir)
{
    this->output_dir = output_dir;
    create_directory(this->output_dir.c_str());
    printf("[config] Output directory = %s\n",this->output_dir.c_str());

    std::string log_filename = this->output_dir + "/output.log";
    this->log_file = fopen(log_filename.c_str(),"w+");
    if (!this->log_file)
    {
        fprintf(stderr,"[config] ERROR! Opening logfile!\n");
        exit(EXIT_FAILURE);
    }
    else
    {
        printf("[config] Log file will be saved at '%s'\n",log_filename.c_str());
    }
}

void Config::print ()
{
    std::string str;
    printf("[config] root point = ( %g, %g, %g )\n",this->root_pos[0],this->root_pos[1],this->root_pos[2]);
    printf("[config] start_diameter = %g\n",this->start_diameter);
    printf("[config] seed = %u\n",this->seed);
    printf("[config] his_offset = %g\n",this->his_offset);
    printf("[config] N_p = %g\n",this->np);
    printf("[config] N_a = %g\n",this->na);
    printf("[config] cost function name = %s\n",this->cost_function_name.c_str());
    str = (this->using_cloud_points) ? "TRUE" : "FALSE";
    printf("[config] using_cloud_points = %s\n",str.c_str());
    str = (this->using_pmj_location) ? "TRUE" : "FALSE";
    printf("[config] using_pmj_location = %s\n",str.c_str());
}
