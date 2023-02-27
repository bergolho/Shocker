// Author: Lucas Berg
// ============================================================================================================
// Program that reads an endocardium mesh with its Local Activation Time (LAT) in VTK format and generate extra
// PMJ points.
// 
// Pre-Requisites:
//  - GCC/G++ 
//  - CMake
//  - Valgrind (for memory leak check ...)
//
//  Features of this program:
//
//   1) Reads a endocardium mesh with its LAT in '.vtk'.
//   2) Generate additional PMJ points using different approaches
//      2.1) Randomly pick a given number of points in the endocardium mesh
//      2.2) Pick a given number of points from the endocardium mesh, but only the ones that are close to a given PMJ set  
//      2.3) Consider a percentage of the endocardium points as the extra PMJs
//   3) Concatenate additional PMJs to a original PMJ set
//    IMPORTANT NOTES: 
//      - The input endocardium mesh must be given in {um} to work (.vtk)
//      - The input PMJ set must be given in {um}. (.vtk)
//      - The output extra PMJs will be given in {um} and {m} with their respective LAT.
//
// How to compile:
//  $ ./recompile_project.sh
// ============================================================================================================


#include "reader/reader.h"

int main (int argc, char *argv[])
{
    if (argc-1 < 1)
    {
        printf("=============================================================================================================================\n");
        printf("Usage:> %s <endocardium_LAT_mesh> [optional_commands]\n",argv[0]);
        printf("-----------------------------------------------------------------------------------------------------------------------------\n");
        printf("<endocardium_LAT_mesh> = Input endocardium mesh in VTU format\n");
        printf("-----------------------------------------------------------------------------------------------------------------------------\n");
        printf("[optional_commands]\n");
        printf("-n <num_pmjs> = Filter the new PMJs by number of points (random)\n");
        printf("-p <percentage> = Filter the new PMJs by percentage of the total number of endocardium points (random)\n");
        printf("-r <original_pmjs_filename> <num_pmjs> = Filter the new PMJs by using a given set of PMJ points. The new PMJs\n");
        printf("                                         will be taken by using a region radius around the given PMJs.\n");
        printf("     * You must supply another input filename with PMJs in VTK format\n");
        printf("     * You must supply a number of PMJs\n");
        printf("=============================================================================================================================\n");
        exit(EXIT_FAILURE);
    }

    int c;
    bool generateByPoints_flag = false;
    bool generateByPercentage_flag = false;
    bool generateByReference_flag = false;
    char *cvalue = NULL;
    opterr = 0;
    while ((c = getopt (argc, argv, "n:p:r:")) != -1)
    {
        switch (c)
        {
            case 'n':
                generateByPoints_flag = true;
                cvalue = optarg;
                break;
            case 'p':
                generateByPercentage_flag = true;
                cvalue = optarg;
                break;
            case 'r':
                generateByReference_flag = true;
                cvalue = optarg;
                break;
            case '?':
                if (optopt == 'n' || optopt == 'p' || optopt == 'r')
                    fprintf (stderr, "Option -%c requires an argument.\n", optopt);
                else if (isprint (optopt))
                    fprintf (stderr, "Unknown option `-%c'.\n", optopt);
                else
                    fprintf (stderr,"Unknown option character `\\x%x'.\n",optopt);
                return 1;
            default:
                abort ();
        }
    }

    std::string endocardium_filename = argv[optind];
    std::vector<Point> endocardium_points;
    printf("[!] Reading endocardium mesh ...\n");
    read_points_from_vtk(endocardium_filename,endocardium_points);

    if (generateByPoints_flag)
    {
        printf("[!] Filtering the PMJs by number ...\n");
        uint32_t num_pmjs = atoi(cvalue);
        std::vector<Point> new_pmjs = filter_points(endocardium_points,num_pmjs);
        write_points_to_vtk("outputs/extra_pmjs_um.vtk",new_pmjs,1);    // um
        write_points_to_vtk("outputs/extra_pmjs_m.vtk",new_pmjs,1e-6);  // m
    }   
    if (generateByPercentage_flag)
    {
        printf("[!] Filtering the PMJs by percentage ...\n");
        double percentage = atof(cvalue);
        uint32_t num_pmjs = endocardium_points.size()*percentage*0.01;
        std::vector<Point> new_pmjs = filter_points(endocardium_points,num_pmjs);
        write_points_to_vtk("outputs/extra_pmjs_um.vtk",new_pmjs,1);    // um
        write_points_to_vtk("outputs/extra_pmjs_m.vtk",new_pmjs,1e-6);  // m
    } 
    if (generateByReference_flag)
    {
        printf("[!] Filtering the PMJs using reference file ...\n");
        std::string original_pmjs_filename = cvalue;
        std::vector<Point> pmj_points;
        read_points_from_vtk(original_pmjs_filename,pmj_points);
        
        uint32_t num_pmjs = atoi(argv[optind+1]);
        std::vector<Point> new_pmjs = filter_points_with_file(endocardium_points,pmj_points,num_pmjs);
        concatenate_points(new_pmjs,pmj_points);
        write_points_to_vtk("outputs/extra_pmjs_um.vtk",new_pmjs,1);    // um
        write_points_to_vtk("outputs/extra_pmjs_m.vtk",new_pmjs,1e-6);  // m
    }

    return 0;
}