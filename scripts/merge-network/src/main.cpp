// Author: Lucas Berg
// =============================================================================================
// Program that reads the 3 Purkinje networks (His-Bundle, LV and RV) and merge them into a single 
//   network.
// 
// Pre-Requisites:
//  - GCC/G++ 
//  - CMake
//  - Valgrind (for memory leak check ...)
//
//  Features of this program:
//
//   1) Read 3 Purkinje networks in '.vkt' format and merge them together into a single '.vtk' file
//   2) Write the merged network in '.vtk' format  
//
// How to compile:
//  $ ./recompile_project.sh
// =============================================================================================

#include "reader/reader.h"
#include "graph/graph.h"

int main (int argc, char *argv[])
{
    if (argc-1 != 5)
    {
        printf("=========================================================================================================================================================================================\n");
        printf("Usage:> %s <input_filename_1> <input_filename_2> <input_filename_3> <root_positions_filename> <output_filename>\n",argv[0]);
        printf("=========================================================================================================================================================================================\n");
        printf("<input_filename_1> = His-bundle\n");
        printf("<input_filename_2> = Left ventricle\n");
        printf("<input_filename_3> = Right ventricle\n");
        printf("<root_positions_filename> = Root positions for the LV and RV trees\n");
        printf("<output_filename> = Merged network\n");
        printf("-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n");
        printf("Example:> %s inputs/reference_his_um.vtk inputs/reference_LV_um.vtk inputs/reference_RV_um.vtk inputs/reference_root_positions.txt outputs/reference_LVRV_um.vtk\n",argv[0]);
        printf("=========================================================================================================================================================================================\n");
	    exit(EXIT_FAILURE);
    }

    std::string his_filename = argv[1];
    VTK_Reader *his_reader = new VTK_Reader(his_filename);
    Graph *his_network = nullptr;
    if (!his_reader->has_point_data)
        his_network = new Graph(his_reader->the_points,his_reader->the_lines);
    else
        his_network = new Graph(his_reader->the_points,his_reader->the_lines,his_reader->the_point_scalars);
    //his_network->print();

    std::string lv_filename = argv[2];
    VTK_Reader *lv_reader = new VTK_Reader(lv_filename);
    Graph *lv_network = nullptr;
    if (!lv_reader->has_point_data)
        lv_network = new Graph(lv_reader->the_points,lv_reader->the_lines);
    else
        lv_network = new Graph(lv_reader->the_points,lv_reader->the_lines,lv_reader->the_point_scalars);
    //lv_network->print();

    std::string rv_filename = argv[3];
    VTK_Reader *rv_reader = new VTK_Reader(rv_filename);
    Graph *rv_network = nullptr;
    if (!rv_reader->has_point_data)
        rv_network = new Graph(rv_reader->the_points,rv_reader->the_lines);
    else
        rv_network = new Graph(rv_reader->the_points,rv_reader->the_lines,rv_reader->the_point_scalars);
    //rv_network->print();

    std::vector<Point> root_pos_arr;
    std::string root_pos_filename = argv[4];
    read_root_positions(root_pos_filename,root_pos_arr);

    Graph *output_network = merge_networks(his_network,lv_network,rv_network,root_pos_arr);
    //output_network->print();

    std::string output_filename = argv[5];
    output_network->write_network(output_filename.c_str());
    //output_network->write_LAT(output_filename.c_str());
    //output_network->write_MonoAlg3D(output_filename.c_str());

    return 0;
}
