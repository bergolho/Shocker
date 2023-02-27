#include "local_optimization_config.h"

LocalOptimizationConfig::LocalOptimizationConfig ()
{
    this->params = new std::map<std::string,std::string>();
}

LocalOptimizationConfig::~LocalOptimizationConfig ()
{
    delete this->params;
}

bool LocalOptimizationConfig::get_parameter_value_from_map (const std::string key, std::string *value)
{
    auto it = this->params->find(key);

    if (it != this->params->end())
    {
        *value = it->second;
        return true;
    }
    else
    {
        fprintf(stderr,"[local_optimization] Not found \"%s\" ... using default value\n",key.c_str());
        return false;
    }
}

void LocalOptimizationConfig::print ()
{
    printf("Local optimization function name = \"%s\"\n",this->function_name.c_str());
    if (!this->params->empty())
    {
        printf("Parameters:\n");
        for (auto it = this->params->begin(); it != this->params->end(); ++it)
            printf("\t%s = %s\n",it->first.c_str(),it->second.c_str());
    }
}