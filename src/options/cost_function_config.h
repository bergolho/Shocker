//
// Created by bergolho on 03/03/19.
//

#ifndef COST_FUNCTION_CONFIG_H
#define COST_FUNCTION_CONFIG_H

#include <cstdio>
#include <cstdlib>
#include <cstdbool>
#include <cstdint>
#include <cstring>
#include <dlfcn.h>

#include <map>
#include <string>

class CostFunctionConfig
{
public:
    std::string function_name;                      // Name of the cost function
    std::map<std::string,std::string> *params;      // Parameters of the cost function
public:
    CostFunctionConfig ();
    ~CostFunctionConfig ();
    bool get_parameter_value_from_map (const std::string key, std::string *value);
    void print ();
};

#endif