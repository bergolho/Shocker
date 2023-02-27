#include "user_options.h"

User_Options::User_Options (const char filename[])
{
    this->use_cloud_points = false;
    this->use_pmj_location = false;
    this->use_initial_network = false;
    this->start_diameter = -1;
    this->np = 400;                             // Default value
    this->na = 400;                             // Default value
    this->seed = 1;                             // Default value
    this->characteristic_length = 5000;         // Default value {um}

    this->cost_function_config = nullptr;
    this->cloud_config = nullptr;
    this->pmj_config = nullptr;

    read_config_file(filename);

    //print();
}

User_Options::~User_Options ()
{
    if (this->cost_function_config)
        delete this->cost_function_config;
    if (this->pmj_config)
        delete this->pmj_config;
    if (this->cloud_config)
        delete this->cloud_config;
}

void User_Options::read_config_file (const char filename[])
{
    printf("%s\n",PRINT_LINE);
    printf("[user_options] Reading configuration file:> \"%s\"\n",filename);

    // Open the config file for reading
    FILE *file = fopen(filename,"r");
    if (!file)
    {
        fprintf(stderr,"[-] Error reading configuration file '%s'\n",filename);
        exit(EXIT_FAILURE);
    }
    
    // Here we parse the config file
    if(ini_parse(filename, parse_config_file, this) < 0) 
    {
        fprintf(stderr, "Error: Can't load the config file %s\n", filename);
        exit(EXIT_FAILURE);
    }

    printf("%s\n",PRINT_LINE);

    fclose(file);

    // DEBUG
    //print();
}

int parse_config_file (void *user, const char *section, const char *name, const char *value)
{
    User_Options *pconfig = (User_Options*)user;

    if (SECTION_STARTS_WITH(MAIN_SECTION))
    {
        if (MATCH_NAME("root_pos"))
        {
            std::string str = value;
            uint32_t first_id = str.find('[');
            uint32_t last_id = str.find(']');
            str = str.substr(first_id+1,last_id-first_id-1);
            
            std::string token;
            std::vector<double> pos;
            std::stringstream ss(str);
            while (std::getline(ss,token,','))
                pos.push_back( std::stof(token) );

            for (uint32_t i = 0; i < 3; i++)
                pconfig->root_pos[i] = pos[i];
        }
        else if (MATCH_NAME("seed"))
        {
            pconfig->seed = (uint32_t)strtol(value, NULL, 10);
        }
        else if (MATCH_NAME("l_d"))
        {
            pconfig->characteristic_length = strtof(value, NULL);
        }
        else if (MATCH_NAME("his_offset"))
        {
            pconfig->his_offset = strtof(value, NULL);
        }
        else if (MATCH_NAME("start_diameter"))
        {
            pconfig->start_diameter = strtof(value, NULL);
        }
        else if (MATCH_NAME("N_p"))
        {
            pconfig->np = (uint32_t)strtol(value, NULL, 10);
        }
        else if (MATCH_NAME("N_a"))
        {
            pconfig->na = (uint32_t)strtol(value, NULL, 10);
        }
        else if (MATCH_NAME("use_initial_network"))
        {
            if (strcmp(value,"true") == 0 || strcmp(value,"yes") == 0)
                pconfig->use_initial_network = true;
            else if (strcmp(value,"false") == 0 || strcmp(value,"no") == 0)
                pconfig->use_initial_network = false;
            else
            {
                fprintf(stderr,"[user_options] Error reading configuration file! Invalid option in \"main\" section\n");
                exit(EXIT_FAILURE);
            }
        }
        else if (MATCH_NAME("initial_network_filename"))
        {
            pconfig->initial_network_filename = value;
        }
    }
    else if (SECTION_STARTS_WITH(SAVE_NETWORK_SECTION))
    {
        if (MATCH_NAME("output_dir"))
        {
            pconfig->output_dir = value;
        }
    }
    else if (SECTION_STARTS_WITH(CLOUD_SECTION))
    {
        if (!pconfig->cloud_config)
        {
            pconfig->cloud_config = new CloudConfig();
        }

        if (MATCH_NAME("use_cloud_points"))
        {
            if (strcmp(value,"true") == 0 || strcmp(value,"yes") == 0)
                pconfig->use_cloud_points = true;
            else if (strcmp(value,"false") == 0 || strcmp(value,"no") == 0)
                pconfig->use_cloud_points = false;
            else
            {
                fprintf(stderr,"[user_options] Error reading configuration file! Invalid option in \"cloud_points\" section\n");
                exit(EXIT_FAILURE);
            }
        }
        else if (MATCH_NAME("cloud_points_filename"))
        {
            pconfig->cloud_config->cloud_filename = value;
        }
        else if (MATCH_NAME("surface_filename"))
        {
            pconfig->cloud_config->surface_filename = value;
        }
        else if (MATCH_NAME("phi"))
        {
            pconfig->cloud_config->phi = strtof(value,NULL);    
        }
        else if (MATCH_NAME("rand_offset"))
        {
            pconfig->cloud_config->rand_offset = (uint32_t)strtol(value, NULL, 10);
        }
    }
    else if (SECTION_STARTS_WITH(PMJ_SECTION))
    {
        if (!pconfig->pmj_config)
        {
            pconfig->pmj_config = new PMJConfig();
        }

        if (MATCH_NAME("use_pmj_location"))
        {
            if (strcmp(value,"true") == 0 || strcmp(value,"yes") == 0)
                pconfig->use_pmj_location = true;
            else if (strcmp(value,"false") == 0 || strcmp(value,"no") == 0)
                pconfig->use_pmj_location = false;
            else
            {
                fprintf(stderr,"[user_options] Error reading configuration file! Invalid option in \"cloud_points\" section\n");
                exit(EXIT_FAILURE);
            }
        }
        else if (MATCH_NAME("pmj_location_filename"))
        {
            pconfig->pmj_config->using_pmj = true;
            pconfig->pmj_config->location_filename = value;
        }
        else if (MATCH_NAME("pmj_connection_rate"))
        {
            pconfig->pmj_config->connection_rate = (uint32_t)strtol(value, NULL, 10);
        }
        else if (MATCH_NAME("lat_error_tolerance"))
        {
            pconfig->pmj_config->lat_error_tolerance = strtof(value, NULL);
        }
    }
    else if (SECTION_STARTS_WITH(COST_FUNCTION_SECTION))
    {
        if (!pconfig->cost_function_config)
        {
            pconfig->cost_function_config = new CostFunctionConfig();
        }

        if (MATCH_NAME("function_name"))
        {
            pconfig->cost_function_config->function_name = value;
        }
        else
        {
            std::string key(name);
            pconfig->cost_function_config->params->insert(std::pair<std::string,std::string>(key,value));
        }
    }
    return 1;
}

void User_Options::print ()
{
    printf("********************* user_options *********************\n");
    printf("seed = %u\n",this->seed);
    printf("root_pos[3] = { %g , %g , %g }\n",this->root_pos[0],this->root_pos[1],this->root_pos[2]);
    printf("start_diameter = %g μm\n",this->start_diameter);
    printf("N_p = %g\n",this->np);
    printf("N_a = %g\n",this->na);
    printf("output_dir = %s\n",this->output_dir.c_str());

    printf("%s\n",PRINT_DOTS);
    if (this->use_cloud_points)
        printf("cloud_points_filename = %s\n",this->cloud_config->cloud_filename.c_str());
    else
        printf("cloud_points_filename = NULL\n");
    if (this->use_initial_network)
        printf("Initial network filename = %s\n",this->initial_network_filename.c_str());
    else
        printf("Initial network filename = NULL\n");
    printf("%s\n",PRINT_DOTS);
    this->cost_function_config->print();
    printf("%s\n",PRINT_DOTS);
    if (this->use_pmj_location)
        this->pmj_config->print();
    else
        printf("PMJ location = FALSE\n");
    printf("********************************************************\n");
}