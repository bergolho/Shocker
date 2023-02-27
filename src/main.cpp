// ========================================================
// Author: Lucas Berg
// ========================================================
// Shocker: A program to generate Purkinje networks.
// ........................................................
// Compile:
//      $ ./recompile_project.sh
// Execute:
//      $ ./bin/Shocker inputs/benchmark_2d.ini
// ========================================================

#include <iostream>
#include <cstdio>
#include <cstdlib>

#include "shocker/shocker.h"

using namespace std;

int main (int argc, char *argv[])
{
    if (!(argc-1 == 1))
    {
        usage(argv[0]);
        exit(EXIT_FAILURE);
    }
    else
    {
        User_Options *options = new User_Options(argv[1]);
        Shocker *shocker = new Shocker(options);

        shocker->grow();
        //shocker->write_full_network_to_vtk();
        
        delete shocker;
        delete options;
    }
    return 0;
}
