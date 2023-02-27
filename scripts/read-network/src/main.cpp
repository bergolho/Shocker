// Author: Lucas Berg
// ============================================================================================================
// Program that reads a Purkinje network stored in .vtk format and outputs some 
// information about it.
// 
// Pre-Requisites:
//  - GCC/G++ 
//  - CMake
//  - Valgrind (for memory leak check ...)
//
//  Features of this program:
//
//   1) Read a Purkinje network in '.vtk' and returns geometric information 
//      (segment length, bifurcation angles and number of segments).
//   2) Generate the LAT map of the given Purkinje network using a constant reference conduction velocity
//   3) Write the positions and the LAT from the terminals of the input Purkinje network.
//   4) Write a percentage of the terminals from the input Purkinje network as active PMJs.
//   5) Compute the errors (RMS, RMSE, eps<2ms, eps<5ms, min.LAT, max.LAT, max.error) between the input
//      Purkinje network and a second one passed as an argument.
//   6) Compute the errors (RMS, RMSE, eps<2ms, eps<5ms, min.LAT, max.LAT, max.error) between the input
//      Purkinje network and a set of PMJs.
//   7) Write the minimum Purkinje network that links the active PMJs from the given input Purkinje network.
//   8) Convert an input Purkinje network to the Shocker initial network topology format.
//   9) Convert an input Purkinje network to a polyline structure
//    IMPORTANT NOTES: 
//      - The input Purkinje network should be given in {um} to work
//      - The output initial network topology file will set all terminal points as ACTIVES
//      - The output initial network topology file will convert all the points to {m}
//      - A default conduction velocity of 1.9m/s will be set to all segments
//      - All segments will have the 'can_write' set to True
//      - The root node of the Purkinje networks should have their index equal to zero
//      - The polyline conversion only works if the input Purkinje network is a single geodesic pathway
//
// How to compile:
//  $ ./recompile_project.sh
// ============================================================================================================

#include "reader/reader.h"
#include "graph/graph.h"

int main (int argc, char *argv[])
{
    if (argc-1 < 2)
    {
        printf("=============================================================================================================================\n");
        printf("Usage:> %s <input_purkinje_network> <input_root_position> [optional_commands]\n",argv[0]);
        printf("-----------------------------------------------------------------------------------------------------------------------------\n");
        printf("<input_purkinje_network> = Input Purkinje network\n");
        printf("<input_root_position> = Input with the root coordinate\n");
        printf("-----------------------------------------------------------------------------------------------------------------------------\n");
        printf("[optional_commands]\n");
        printf("-g = Compute geometrics for the input Purkinje network\n");
        printf("-l = Write the LAT map for the input Purkinje network\n");
        printf("-t = Write the terminal locations of the input Purkinje network\n");
        printf("-c = Convert the input Purkinje network to the Shocker initial network topology format\n");
        printf("-i = Write graph information about the input Purkinje network\n");
        printf("-p = Write a polyline of a given input Purkinje network\n");
        printf("-a <percentage> = Write a percentage of the terminals from the input Purkinje network as active PMJs\n");
        printf("-m <pmj_filename> = Write the minimum network that links the active PMJs from the input network\n");
        printf("-e <input_filename> = Compute the LAT errors between the input Purkinje network and another one\n");
        printf("     You must supply another input filename with a Purkinje network in the same format\n");
        printf("-r <input_filename> = Compute the LAT errors between the input Purkinje network and a set of PMJ points\n");
        printf("     You must supply another input filename with the coordinates and LAT of the PMJs\n");
        printf("=============================================================================================================================\n");
        exit(EXIT_FAILURE);
    }

    int c;
    bool compGeometrics_flag = false;
    bool compErrors_flag = false;
    bool compErrorsWithPMJ_flag = false;
    bool convertToShocker_flag = false;
    bool writeLAT_flag = false;
    bool writeTerminals_flag = false;
    bool writeActivePMJ_flag = false;
    bool writeMinimumTree_flag = false;
    bool writeGraphInfo_flag = false;
    bool writePolyne_flag = false;
    char *cvalue = NULL;
    opterr = 0;
    while ((c = getopt (argc, argv, "gtlcipe:a:m:r:")) != -1)
    {
        switch (c)
        {
            case 'g':
                compGeometrics_flag = true;
                break;
            case 'l':
                writeLAT_flag = true;
                break;
            case 't':
                writeTerminals_flag = true;
                break;
            case 'i':
                writeGraphInfo_flag = true;
                break;
            case 'c':
                convertToShocker_flag = true;
                break;
            case 'p':
                writePolyne_flag = true;
                break;
            case 'a':
                writeActivePMJ_flag = true;
                cvalue = optarg;
                break;
            case 'm':
                writeMinimumTree_flag = true;
                cvalue = optarg;
                break;
            case 'e':
                compErrors_flag = true;
                cvalue = optarg;
                break;
            case 'r':
                compErrorsWithPMJ_flag = true;
                cvalue = optarg;
                break;
            case '?':
                if (optopt == 'e' || optopt == 'a' || optopt == 'm')
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

    // Reads and builds the input Purkinje network
    std::string pk_filename = argv[optind];
    std::string root_filename = argv[optind+1];
    VTK_Reader *reader = new VTK_Reader(pk_filename);
    Graph *network = new Graph(reader->the_points,reader->the_lines,reader->the_celldata,root_filename);

    if (compGeometrics_flag)
    {
        printf("[+] Computing geometric values ...\n");
        network->compute_geometrics();
    }
    if (writeTerminals_flag)
    {
        printf("[+] Writing terminal locations ...\n");
        network->write_terminals("outputs/terminals.vtk");
    }
    if (writeLAT_flag)
    {
        printf("[+] Writing LAT map ...\n");
        network->write_LAT("outputs/network.vtk");
    }
    if (writeActivePMJ_flag)
    {
        double percentage = atof(cvalue);
        printf("[+] Writing active PMJs ...\n");
        network->write_active_pmjs("outputs/pmjs.vtk",percentage);
    }
    if (writeMinimumTree_flag)
    {
        std::string pmj_filename = cvalue;
        printf("[+] Writing minimum network ...\n");
        network->write_minimum_network(pmj_filename,"outputs/minimum_network.vtk");
    }
    if (writeGraphInfo_flag)
    {
        printf("[+] Writing graph information ...\n");
        network->write_graph_info("outputs/graph_information.vtk");
    }
    if (convertToShocker_flag)
    {
        printf("[+] Converting to Shocker initial network topology ...\n");
        network->convert_to_initial_network_topology("outputs/network_topology.txt");
    }
    if (writePolyne_flag)
    {
        printf("[+] Writing polyline ...\n");
        network->write_polyline("outputs/polyline.vtk");
    }
    if (compErrors_flag)
    {
        std::string ref_filename = cvalue;
        VTK_Reader *ref_reader = new VTK_Reader(ref_filename);
        Graph *ref_network = new Graph(ref_reader->the_points,ref_reader->the_lines,ref_reader->the_celldata,root_filename);
        
        network->compute_error(ref_network);
    }
    if (compErrorsWithPMJ_flag)
    {
        std::string pmj_filename = cvalue;
        std::vector<PMJ> pmjs;
        read_pmjs(pmj_filename,pmjs);

        network->compute_error(pmjs);
    }
    
    return 0;
}
