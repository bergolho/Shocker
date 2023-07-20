// =============================================================================
// Author: Lucas Berg
// Last changed: 19/07/2023
// =============================================================================
// Shocker: A program to generate patient-specific Purkinje networks.
// ........................................................,,,,,,,,,,,,,,,,,,,,,
// Compile:
//      $ ./recompile_project.sh
// Execute:
//      $ ./bin/Shocker inputs/simplified_LV.ini
//   or
//      $ ./bin/Shocker inputs/simplified_RV.ini
// ........................................................,,,,,,,,,,,,,,,,,,,,,
// To merge the generated LV and RV into an LVRV the user needs to manually
// construct a His-Bundle in .vtk format. There is an example for the Simplified
// mesh inside the "/scripts/merge-network" folder.
// =============================================================================

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
