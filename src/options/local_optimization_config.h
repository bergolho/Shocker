#ifndef LOCAL_OPTIMIZATION_CONFIG_H
#define LOCAL_OPTIMIZATION_CONFIG_H

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>

#include <string>
#include <map>

class LocalOptimizationConfig
{
public:
    std::string function_name;                      // Local optimization function name
    std::map<std::string,std::string> *params;      // Parameters of the local optimization function
public:
    LocalOptimizationConfig ();
    ~LocalOptimizationConfig ();
    bool get_parameter_value_from_map (const std::string key, std::string *value);
    void print ();
};

// Auxiliary functions

#endif