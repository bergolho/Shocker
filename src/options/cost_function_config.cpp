#include "cost_function_config.h"

CostFunctionConfig::CostFunctionConfig ()
{
    this->params = new std::map<std::string,std::string>();
}

CostFunctionConfig::~CostFunctionConfig ()
{
    delete this->params;
}

bool CostFunctionConfig::get_parameter_value_from_map (const std::string key, std::string *value)
{
    auto it = this->params->find(key);

    if (it != this->params->end())
    {
        *value = it->second;
        return true;
    }
    else
    {
        fprintf(stderr,"Not found \"%s\" ... using default value\n",key.c_str());
        return false;
    }
}

void CostFunctionConfig::print ()
{
    printf("Cost function name = \"%s\"\n",this->function_name.c_str());
    if (!this->params->empty())
    {
        printf("Parameters:\n");
        for (auto it = this->params->begin(); it != this->params->end(); ++it)
            printf("\t%s = %s\n",it->first.c_str(),it->second.c_str());
    }
}