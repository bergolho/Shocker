// Author: Lucas Berg
// =============================================================================================
// Program that reads 2 Purkinje networks and concatenates them into a single '.vtk' file.
// 
// Pre-Requisites:
//  - GCC/G++ 
//  - CMake
//  - Valgrind (for memory leak check ...)
//
//  Features of this program:
//
//   1) Read 2 Purkinje networks in '.vkt' format and concatenate the points and lines into a 
//      single '.vtk' file
//   2) Write the concatenated network in '.vtk' format  
//
// How to compile:
//  $ ./recompile_project.sh
// =============================================================================================

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <cmath>
#include <vector>

uint32_t total_number_nodes;
uint32_t total_number_edges;

struct node
{
    uint32_t id;
    double x, y, z;
    std::vector< std::pair<uint32_t,double> > edges;
};

void new_purkinje_network(const char filename[], std::vector< struct node > &pk)
{
    FILE *file = fopen(filename,"r");

    char str[200];
    while (fscanf(file,"%s",str) != EOF)
    {
        if (strcmp(str,"POINTS") == 0) break;
    }
    int num_nodes;
    fscanf(file,"%d %s",&num_nodes,str);

    for (int i = 0; i < num_nodes; i++)
    {
        double x, y, z;
        fscanf(file,"%lf %lf %lf",&x,&y,&z);

        struct node n;
        n.id = i;
        n.x = x;
        n.y = y;
        n.z = z;
        pk.push_back(n);
    }

    while (fscanf(file,"%s",str) != EOF)
    {
        if (strcmp(str,"LINES") == 0) break;
    }    
    int num_edges, trash;
    std::vector< std::pair<int,int> > lines;
    fscanf(file,"%d %d",&num_edges,&trash);
    for (int i = 0; i < num_edges; i++)
    {
        int src, dest;
        fscanf(file,"%d %d %d",&trash,&src,&dest);

        pk[src].edges.push_back( std::make_pair(dest,0) );
        lines.push_back( std::make_pair(src,dest) );
    }
    while (fscanf(file,"%s",str) != EOF)
    {
        if (strcmp(str,"default") == 0) break;
    }
    
    for (int i = 0; i < num_edges; i++)
    {
        double value;
        fscanf(file,"%lf",&value);

        int src = lines[i].first;
        int dest = lines[i].second;
        for (uint32_t j = 0; j < pk[src].edges.size(); j++)
        {
            if (pk[src].edges[j].first == dest)
                pk[src].edges[j].second = value;
        }
    }

    fclose(file);
} 

void print_network (std::vector< struct node > pk)
{
    for (uint32_t i = 0; i < pk.size(); i++)
    {
        printf("|| %d %g %g %g || ",pk[i].id,pk[i].x,pk[i].y,pk[i].z);
        for (uint32_t j = 0; j < pk[i].edges.size(); j++)
        {
            printf(" --> || %d ||",pk[i].edges[j].first);
        }
        printf("\n");
    }
}

uint32_t get_root_index (std::vector< struct node > lv_network, const double root[3])
{
    double min_dist = __DBL_MAX__;
    uint32_t min_index = 0;

    for (uint32_t i = 0; i < lv_network.size(); i++)
    {
        double p1[3];
        p1[0] = lv_network[i].x;
        p1[1] = lv_network[i].y;
        p1[2] = lv_network[i].z;

        double dist = sqrt(pow(p1[0]-root[0],2) + pow(p1[1]-root[1],2) + pow(p1[2]-root[2],2));
        if (dist < min_dist)
        {
            min_dist = dist;
            min_index = i;
        }
    }
    return min_index;
}

void merge_purkinje_networks (std::vector< struct node > network_1, std::vector< struct node > network_2,\
                            std::vector< struct node > &output_network)
{
    uint32_t offset = 0;
    for (uint32_t i = 0; i < network_1.size(); i++)
    {
        struct node n;
        n.id = i;
        n.x = network_1[i].x;
        n.y = network_1[i].y;
        n.z = network_1[i].z;

        output_network.push_back(n);
    }
    offset = network_1.size();
    for (uint32_t i = 0; i < network_2.size(); i++)
    {
        struct node n;
        n.id = i + offset;
        n.x = network_2[i].x;
        n.y = network_2[i].y;
        n.z = network_2[i].z;

        output_network.push_back(n);
    }
    total_number_nodes = network_1.size() + network_2.size();

    // Edges for the first network
    offset = 0;
    total_number_edges = 0;
    for (uint32_t i = 0; i < network_1.size(); i++)
    {
        int src = i;
        for (uint32_t j = 0; j < network_1[i].edges.size(); j++)
        {
            uint32_t dest = network_1[i].edges[j].first;
            double radius = network_1[i].edges[j].second;
            output_network[src].edges.push_back( std::make_pair(dest,radius) );
            total_number_edges++;
        }
    }

    // Set edges for the LV
    offset = network_1.size();
    for (uint32_t i = 0; i < network_2.size(); i++)
    {
        uint32_t src = i;
        for (uint32_t j = 0; j < network_2[i].edges.size(); j++)
        {
            uint32_t dest = network_2[i].edges[j].first;
            double radius = network_2[i].edges[j].second;
            output_network[src+offset].edges.push_back( std::make_pair(dest+offset,radius) );
            total_number_edges++;
        }
    }
}

void write_network_to_vtk(const char filename[], std::vector< struct node > new_network)
{
    FILE *file = fopen(filename,"w+");

    fprintf(file,"# vtk DataFile Version 3.0\n");
    fprintf(file,"Tree\n");
    fprintf(file,"ASCII\n");
    fprintf(file,"DATASET POLYDATA\n");
    fprintf(file,"POINTS %lu float\n",new_network.size());
    for (uint32_t i = 0; i < new_network.size(); i++)
        fprintf(file,"%g %g %g\n",new_network[i].x,new_network[i].y,new_network[i].z);
    fprintf(file,"LINES %u %u\n",total_number_edges,total_number_edges*3);
    for (uint32_t i = 0; i < new_network.size(); i++)
    {
        for (uint32_t j = 0; j < new_network[i].edges.size(); j++)
        {
            fprintf(file,"2 %u %u\n",i,new_network[i].edges[j].first);
        }
    }
    
    fprintf(file,"CELL_DATA %u\n",total_number_edges);
    fprintf(file,"SCALARS radius float\n");
    fprintf(file,"LOOKUP_TABLE default\n");
    for (uint32_t i = 0; i < new_network.size(); i++)
    {
        for (uint32_t j = 0; j < new_network[i].edges.size(); j++)
        {
            fprintf(file,"%g\n",new_network[i].edges[j].second);
        }
    }

    //fprintf(file,"VERTICES %u %u\n",new_network.size(),new_network.size()*2);
    //for (uint32_t i = 0; i < new_network.size(); i++)
    //    fprintf(file,"1 %u\n",i);

    fclose(file);
}

void write_network_terminals_to_vtk (std::vector< struct node > network)
{
    std::vector<uint32_t> term_ids;
    for (uint32_t i = 0; i < network.size(); i++)
    {
        if (network[i].edges.size() == 0)
            term_ids.push_back(i);
    }

    FILE *file = fopen("outputs/network_terminals.vtk","w+");
    fprintf(file,"# vtk DataFile Version 3.0\n");
    fprintf(file,"Tree\n");
    fprintf(file,"ASCII\n");
    fprintf(file,"DATASET POLYDATA\n");
    fprintf(file,"POINTS %lu float\n",term_ids.size());
    for (uint32_t i = 0; i < term_ids.size(); i++)
    {
        uint32_t id = term_ids[i];
        fprintf(file,"%g %g %g\n",network[id].x,network[id].y,network[id].z);
    }
    fprintf(file,"VERTICES %u %u\n",term_ids.size(),term_ids.size()*2);
    for (uint32_t i = 0; i < term_ids.size(); i++)
        fprintf(file,"1 %u\n",i);
    fclose(file);
}

int main (int argc, char *argv[])
{
    if (argc-1 != 3)
    {
        printf("========================================================================================================\n");
        printf("Usage:> %s <input_filename_1> <input_filename_2> <output_filename>\n",argv[0]);
        printf("========================================================================================================\n");
        printf("<input_filename_1> = Purkinje network 1\n");
        printf("<input_filename_2> = Purkinje network 2\n");
        printf("<output_filename> = Concatenated network\n");
        printf("========================================================================================================\n");
	    exit(EXIT_FAILURE);
    }

    char *network_filename_1 = argv[1];
    char *network_filename_2 = argv[2];
    char *output_network_filename = argv[3];

    std::vector< struct node > network_1;
    new_purkinje_network(network_filename_1,network_1);
    //print_network(network_1);

    //write_network_terminals_to_vtk(network_1);

    std::vector< struct node > network_2;
    new_purkinje_network(network_filename_2,network_2);
    //print_network(network_2);

    std::vector< struct node > output_network;
    merge_purkinje_networks(network_1,network_2,output_network);
    //print_network(output_network);

    write_network_to_vtk(output_network_filename,output_network);

    return 0;
}
