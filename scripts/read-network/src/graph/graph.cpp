#include "graph.h"

int counter = 0;

Graph::Graph ()
{

}

Graph::Graph (std::vector<Point> points, std::vector<Line> lines, std::vector<std::vector<double>> cell_data, std::string root_filename)
{
    double root_pos[3];
    read_root_position(root_filename,root_pos);

    uint32_t num_nodes = points.size();
    uint32_t num_edges = lines.size();

    this->list_nodes.assign(num_nodes,Node());
    for (uint32_t i = 0; i < num_nodes; i++)
    {
        uint32_t id = points[i].id;
        double x = points[i].x;
        double y = points[i].y;
        double z = points[i].z;

        this->list_nodes[i].id = id;
        this->list_nodes[i].x = x;
        this->list_nodes[i].y = y;
        this->list_nodes[i].z = z;
    }


    for (uint32_t i = 0; i < num_edges; i++)
    {
        uint32_t src = lines[i].src;
        uint32_t dest = lines[i].dest;
        double length = calc_norm(points[src].x,points[src].y,points[src].z,\
                                points[dest].x,points[dest].y,points[dest].z);
        
        Edge e(dest,length);
        if (cell_data.size() > 0) {

            double _lat = cell_data[0][i];
            double _diameter = cell_data[1][i];
            double _length = cell_data[2][i];
            double _min_tree = cell_data[3][i];

            e._lat = _lat;
            e._diameter = _diameter;
            e._length = _length;
            e._min_tree = _min_tree;

        }
        this->list_nodes[src].list_edges.push_back(e);
    }

    total_nodes = num_nodes;
    total_edges = num_edges;

    // Run a Dijkstra shortest path to fill the 'dist' and 'parent' arrays
    this->root_index = get_root_position_index(root_pos);
    dijkstra(this->root_index);

    // Compute the LAT of all cells using the reference conduction velocity
    this->min_LAT = __DBL_MAX__;
    compute_activation_times(REF_CV);
    //compute_activation_times_with_diameter();

    // Search for the terminal points and fill the 'terminal_indexes' array
    fill_terminal_indexes(this->root_index);
}

void Graph::depth_first_search_default (const uint32_t src_id)
{
    std::vector<bool> dfs_visited;
    dfs_visited.assign(this->total_nodes,false);

    dfs_default(this->list_nodes[src_id],dfs_visited);
}

void Graph::depth_first_search_topology (const uint32_t src_id)
{
    this->dfs_visited.assign(this->total_nodes,false);
    this->parent.assign(this->total_nodes,-1);
    this->dfs_counter.assign(this->total_nodes,-1);

    Node *src_node = &this->list_nodes[src_id];
    dfs_topology(src_node);
}

void Graph::depth_first_search_branches (const uint32_t src_id, std::vector<double> &the_branches)
{
    this->dfs_visited.assign(this->total_nodes,false);

    uint32_t total_branches = 0;
    double segment_size = 0.0;
    uint32_t flag = 0;

    dfs_branches(this->list_nodes[src_id],the_branches);
    printf("Total branches = %u\n",total_branches);
}

void Graph::dfs_default (Node u, std::vector<bool> &dfs_visited)
{
    uint32_t u_index = u.id;
    dfs_visited[u_index] = true;

    uint32_t num_edges = u.list_edges.size();
    for (uint32_t j = 0; j < num_edges; j++)
    {
        uint32_t v_index = u.list_edges[j].dest_id;
        double length = u.list_edges[j].length;

        if (!dfs_visited[v_index])
            dfs_default(this->list_nodes[v_index],dfs_visited);
    }
}

void Graph::dfs_branches (Node u, std::vector<double> &segments)
{
    uint32_t u_index = u.id;
    uint32_t num_edges = u.list_edges.size();
    this->dfs_visited[u_index] = true;
    if (num_edges == 2)
    {
        
    }

    
    for (uint32_t j = 0; j < num_edges; j++)
    {
        uint32_t v_index = u.list_edges[j].dest_id;
        double length = u.list_edges[j].length;

        if (!this->dfs_visited[v_index])
        {
            dfs_branches(this->list_nodes[v_index],segments);
        }
    }
}

void Graph::dfs_topology (Node *current)
{
    uint32_t u_index = current->id;
    this->dfs_visited[u_index] = true;
    this->dfs_counter[u_index] = counter; counter++;
    
    uint32_t num_edges = current->list_edges.size();
    for (uint32_t j = 0; j < num_edges; j++)
    {
        uint32_t v_index = current->list_edges[j].dest_id;
        double length = current->list_edges[j].length;

        if (!this->dfs_visited[v_index])
        {
            this->parent[v_index] = u_index;
            dfs_topology(&this->list_nodes[v_index]);
        }
    }
}

void Graph::print ()
{
    for (uint32_t i = 0; i < this->total_nodes; i++)
    {
        this->list_nodes[i].print();
    }
}

uint32_t Graph::get_closest_terminal_point (Point p, const bool is_reference)
{
    uint32_t biff_value = (is_reference) ? 1 : 0;
    uint32_t closest_index = 0;
    double closest_dist = __DBL_MAX__;
    for (uint32_t i = 0; i < this->list_nodes.size(); i++)
    {
        //if (this->list_nodes[i].list_edges.size() == biff_value)
        //{
            double dist = calc_norm(this->list_nodes[i].x,this->list_nodes[i].y,this->list_nodes[i].z,\
                                p.x,p.y,p.z);

            if (dist < closest_dist)
            {
                closest_dist = dist;
                closest_index = i;
            }
        //}
    }
    return closest_index;
}

// Using C++/STD library
void Graph::dijkstra (const uint32_t src_id)
{
    uint32_t np = this->total_nodes;
    dist.assign(np,__DBL_MAX__);
    parent.assign(np,-1);
    dist[src_id] = 0.0;

    std::priority_queue< std::pair<double,uint32_t>, std::vector< std::pair<double,uint32_t> >, std::greater< std::pair<double,uint32_t> > > pq;
    pq.push(std::make_pair(0.0,src_id));

    while (!pq.empty())
    {
        std::pair<double,uint32_t> front = pq.top(); pq.pop();
        double d = front.first;
        uint32_t u = front.second;
        if (d > dist[u]) 
            continue;
        
        for (uint32_t i = 0; i < this->list_nodes[u].list_edges.size(); i++)
        {
            uint32_t v = this->list_nodes[u].list_edges[i].dest_id;
            double w = this->list_nodes[u].list_edges[i].length;

            if (dist[u] + w < dist[v])
            {
                dist[v] = dist[u] + w;
                parent[v] = u;
                pq.push(std::make_pair(dist[v],v));
            }
        }
    }

    //for (uint32_t i = 0; i < dist.size(); i++)
    //    printf("Node %u -- Parent = %d -- Dist = %g\n",i,parent[i],dist[i]);
}

// Using PQUEUE library
void Graph::dijkstra_2 (const uint32_t src_id)
{
    // Initialize the shortest distance array
    uint32_t np = this->total_nodes;
    dist.assign(np,__DBL_MAX__);
    dist[src_id] = 0.0;

    pqueue_t *pq;
	node_t   *ns;
	node_t   *n;

    // Initialize the priority queue
	ns = (struct node_t*)malloc(np * sizeof(node_t));
	pq = pqueue_init(np, cmp_pri, get_pri, set_pri, get_pos, set_pos);
	if (!(ns && pq))
    {
        fprintf(stderr,"ERROR! Bad alloc!\n");
        exit(EXIT_FAILURE); 
    }
    ns[src_id].pri = 0.0; ns[src_id].val = src_id; pqueue_insert(pq, &ns[src_id]);

    while ((n = (node_t*)pqueue_pop(pq)))
    {
        double d = n->pri;
        uint32_t u = n->val;
        if (d > dist[u]) 
            continue;
        
        for (uint32_t i = 0; i < this->list_nodes[u].list_edges.size(); i++)
        {
            uint32_t v = this->list_nodes[u].list_edges[i].dest_id;
            double w = this->list_nodes[u].list_edges[i].length;

            if (dist[u] + w < dist[v])
            {
                dist[v] = dist[u] + w;
                ns[v].pri = dist[v]; ns[v].val = v; pqueue_insert(pq, &ns[v]);
            }
        }
    }

	pqueue_free(pq);
	free(ns);

    //for (uint32_t i = 0; i < dist.size(); i++)
    //    printf("Node %u -- Dist = %g\n",i,dist[i]);
}

void Graph::fill_terminal_indexes (const uint32_t root_index)
{
    Node src = this->list_nodes[root_index];
    uint32_t num_nodes = this->list_nodes.size();
    std::vector<bool> dfs_visited;
    dfs_visited.assign(num_nodes,false);

    // Visit only the nodes that are part of the network
    dfs_default(src,dfs_visited);
    
    // Pass through all nodes, but only consider the visited ones
    for (uint32_t i = 0; i < num_nodes; i++)
    {
        Node cur_node = this->list_nodes[i];

        // Check if the node is a terminal and was visited
        if (cur_node.is_terminal() && dfs_visited[i])
            this->terminals_indexes.push_back(i);
    }

    // Sort the terminal indexes by (x,y,z)
    std::vector<Point> sorted_terminals;
    for (uint32_t i = 0; i < this->terminals_indexes.size(); i++)
    {
        uint32_t id = this->terminals_indexes[i];
        Point p(id,this->list_nodes[id].x,this->list_nodes[id].y,this->list_nodes[id].z);
        sorted_terminals.push_back(p);
    }
    std::sort(sorted_terminals.begin(),sorted_terminals.end(),compareXYZ);

    // Update the terminal indexes with the sorted indexes by (x,y,z)
    for (uint32_t i = 0; i < this->terminals_indexes.size(); i++)
        this->terminals_indexes[i] = sorted_terminals[i].id;
}

void Graph::update_reference_closest_terminal_indexes (Graph *ref_network)
{
    uint32_t num_terminals = this->terminals_indexes.size();
    
    for (uint32_t i = 0; i < num_terminals; i++)
    {
        uint32_t closest_id = 0;
        double closest_dist = __DBL_MAX__;

        uint32_t cur_id = this->terminals_indexes[i];
        Node cur_term = this->list_nodes[cur_id];

        for (uint32_t j = 0; j < ref_network->list_nodes.size(); j++)
        {
            Node n = ref_network->list_nodes[j];
            if (n.is_terminal())
            {
                double dist = calc_norm(n.x,n.y,n.z,cur_term.x,cur_term.y,cur_term.z);
                if (dist < closest_dist)
                {
                    closest_dist = dist;
                    closest_id = j;
                }
            }
        }

        ref_network->terminals_indexes[i] = closest_id;
    }

}

void Graph::compute_error (Graph *ref_network)
{
    printf("%u %u\n",this->terminals_indexes.size(),ref_network->terminals_indexes.size());
    assert( this->terminals_indexes.size() == ref_network->terminals_indexes.size() );

    update_reference_closest_terminal_indexes(ref_network);

    // Minimum and maximum LAT
    double min_lat, max_lat;
    double ref_min_lat, ref_max_lat;
    this->compute_min_max_lat(min_lat,max_lat);
    ref_network->compute_min_max_lat(ref_min_lat,ref_max_lat);

    // Maximum error
    double max_error;
    this->compute_max_error(max_error,ref_network);

    // RMSE and RRMSE
    double rmse, rrmse;
    this->compute_rmse_rrmse(ref_network,rmse,rrmse);

    double epsilon_2ms, epsilon_5ms;
    this->compute_epsilon_percentage(ref_network,2,epsilon_2ms);
    this->compute_epsilon_percentage(ref_network,5,epsilon_5ms);

    printf("[INFO] Minimum LAT = %g ms\n",min_lat);
    printf("[INFO] Maximum LAT = %g ms\n",max_lat);
    printf("[INFO] Minimum LAT (reference) = %g ms\n",ref_min_lat);
    printf("[INFO] Maximum LAT (reference) = %g ms\n",ref_max_lat);
    printf("[INFO] Maximum error = %g ms\n",max_error);
    printf("[INFO] RMSE = %g ms\n",rmse);
    printf("[INFO] RRMSE = %g \%\n",rrmse*100);
    printf("[INFO] Epsilon 2ms = %g \%\n",epsilon_2ms);
    printf("[INFO] Epsilon 5ms = %g \%\n",epsilon_5ms);

    write_electric_info_to_file("outputs/electric_info.txt",min_lat,max_lat,max_error,rmse,rrmse,epsilon_2ms,epsilon_5ms);
}

void Graph::compute_error (std::vector<PMJ> pmjs)
{
    printf("%u %u\n",this->terminals_indexes.size(),pmjs.size());
    assert( this->terminals_indexes.size() == pmjs.size() );

    // Minimum and maximum LAT
    double min_lat, max_lat;
    double ref_min_lat, ref_max_lat;
    this->compute_min_max_lat(min_lat,max_lat);
    compute_min_max_lat_pmjs(pmjs,ref_min_lat,ref_max_lat);

    // Maximum error
    double max_error;
    this->compute_max_error(max_error,pmjs);

    // RMSE and RRMSE
    double rmse, rrmse;
    this->compute_rmse_rrmse(pmjs,rmse,rrmse);

    double epsilon_2ms, epsilon_5ms;
    this->compute_epsilon_percentage(pmjs,2,epsilon_2ms);
    this->compute_epsilon_percentage(pmjs,5,epsilon_5ms);

    printf("[INFO] Minimum LAT = %g ms\n",min_lat);
    printf("[INFO] Maximum LAT = %g ms\n",max_lat);
    printf("[INFO] Minimum LAT (reference) = %g ms\n",ref_min_lat);
    printf("[INFO] Maximum LAT (reference) = %g ms\n",ref_max_lat);
    printf("[INFO] Maximum error = %g ms\n",max_error);
    printf("[INFO] RMSE = %g ms\n",rmse);
    printf("[INFO] RRMSE = %g \%\n",rrmse*100);
    printf("[INFO] Epsilon 2ms = %g \%\n",epsilon_2ms);
    printf("[INFO] Epsilon 5ms = %g \%\n",epsilon_5ms);
//
    //write_electric_info_to_file("outputs/electric_info.txt",min_lat,max_lat,max_error,rmse,rrmse,epsilon_2ms,epsilon_5ms);
}

void Graph::compute_activation_times (const double cv)
{
    this->lat.clear();
    for (uint32_t i = 0; i < this->total_nodes; i++)
    {
        double lat = (this->dist[i] != __DBL_MAX__) ? this->dist[i]/cv : -1;
        this->lat.push_back(lat);
    }
}

void Graph::compute_min_max_lat (double &min_value, double &max_value)
{
    uint32_t num_terminals = this->terminals_indexes.size();
    max_value = __DBL_MIN__;
    min_value = __DBL_MAX__;
    for (uint32_t i = 0; i < num_terminals; i++)
    {
        uint32_t id = this->terminals_indexes[i];
        if (this->lat[id] > max_value)
        {
            max_value = this->lat[id];
            //max_id = index;
        }
        if (this->lat[id] < min_value)
        {
            min_value = this->lat[id];
            //min_id = index;
        }
    }
}

void Graph::compute_max_error (double &max_error, Graph *ref_network)
{
    uint32_t num_terminals = this->terminals_indexes.size();
    max_error = __DBL_MIN__;
    for (uint32_t i = 0; i < num_terminals; i++)
    {
        uint32_t id = this->terminals_indexes[i];
        uint32_t ref_id = ref_network->terminals_indexes[i];
        double error = fabs(this->lat[id] - ref_network->lat[ref_id]);
        if (error > max_error)
            max_error = error;
    }
}

void Graph::compute_max_error (double &max_error, std::vector<PMJ> ref_pmjs)
{
    uint32_t num_terminals = this->terminals_indexes.size();
    max_error = __DBL_MIN__;
    for (uint32_t i = 0; i < num_terminals; i++)
    {
        uint32_t id = this->terminals_indexes[i];
        uint32_t closest_id = get_closest_pmj(this->list_nodes[id].x,this->list_nodes[id].y,this->list_nodes[id].z,ref_pmjs);
        
        double error = fabs(this->lat[id] - ref_pmjs[closest_id].lat);
        if (error > max_error)
            max_error = error;
    }
}

void Graph::compute_rmse_rrmse (Graph *ref_network, double &rmse, double &rrmse)
{
    uint32_t num_terminals = this->terminals_indexes.size();
    double sum_num = 0.0;
    double sum_den = 0.0;
    for (uint32_t i = 0; i < num_terminals; i++)
    {
        uint32_t id = this->terminals_indexes[i];
        uint32_t ref_id = ref_network->terminals_indexes[i];
        double error = fabs(this->lat[id] - ref_network->lat[ref_id]);

        sum_num += powf(error,2);
        sum_den += powf(ref_network->lat[ref_id],2);
    }    
    double l2_norm = sqrt(sum_den);
    rmse = sqrt(sum_num/(double)num_terminals);
    rrmse = sqrt(sum_num/sum_den);
}

void Graph::compute_rmse_rrmse (std::vector<PMJ> ref_pmjs, double &rmse, double &rrmse)
{
    uint32_t num_terminals = this->terminals_indexes.size();
    double sum_num = 0.0;
    double sum_den = 0.0;
    for (uint32_t i = 0; i < num_terminals; i++)
    {
        uint32_t id = this->terminals_indexes[i];
        uint32_t closest_id = get_closest_pmj(this->list_nodes[id].x,this->list_nodes[id].y,this->list_nodes[id].z,ref_pmjs);

        double error = fabs(this->lat[id] - ref_pmjs[closest_id].lat);

        sum_num += powf(error,2);
        sum_den += powf(ref_pmjs[closest_id].lat,2);
    }    
    double l2_norm = sqrt(sum_den);
    rmse = sqrt(sum_num/(double)num_terminals);
    rrmse = sqrt(sum_num/sum_den);
}

void Graph::compute_epsilon_percentage (Graph *ref_network, const double epsilon, double &percentage)
{
    uint32_t counter = 0;
    uint32_t num_terminals = this->terminals_indexes.size();
    for (uint32_t i = 0; i < num_terminals; i++)
    {
        uint32_t id = this->terminals_indexes[i];
        uint32_t ref_id = ref_network->terminals_indexes[i];
        double error = fabs(this->lat[id] - ref_network->lat[ref_id]);
        
        if (error < epsilon) counter++;
    }    
    percentage = (double)counter / (double)num_terminals * 100.0;
}

void Graph::compute_epsilon_percentage (std::vector<PMJ> ref_pmjs, const double epsilon, double &percentage)
{
    uint32_t counter = 0;
    uint32_t num_terminals = this->terminals_indexes.size();
    for (uint32_t i = 0; i < num_terminals; i++)
    {
        uint32_t id = this->terminals_indexes[i];
        uint32_t closest_id = get_closest_pmj(this->list_nodes[id].x,this->list_nodes[id].y,this->list_nodes[id].z,ref_pmjs);
        
        double error = fabs(this->lat[id] - ref_pmjs[closest_id].lat);
        
        if (error < epsilon) counter++;
    }    
    percentage = (double)counter / (double)num_terminals * 100.0;
}

void Graph::compute_geometrics ()
{
    std::vector<double> segments;
    get_segment_length(segments);
    write_data_to_file("outputs/all_segment_length.dat",segments);

    std::vector<double> branches;
    get_branches(branches);
    write_data_to_file("outputs/all_branches.dat",branches);

    std::vector<double> angles;
    get_bifurcation_angles(angles);
    write_data_to_file("outputs/all_bifurcation_angles.dat",angles);

    double mean_segment_length, std_segment_length;
    compute_mean_std(segments,mean_segment_length,std_segment_length);

    double mean_branch_length, std_branch_length;
    compute_mean_std(branches,mean_branch_length,std_branch_length);

    double mean_biff_angle, std_biff_angle;
    compute_mean_std(angles,mean_biff_angle,std_biff_angle);

    printf("[INFO] Total number of segment = %u\n",segments.size());
    printf("[INFO] Segment length = %.2lf +/- %.2lf mm\n",mean_segment_length*UM_TO_MM,std_segment_length*UM_TO_MM);
    printf("[INFO] Total number of branches = %u\n",branches.size());
    printf("[INFO] Branch length = %.2lf +/- %.2lf mm\n",mean_branch_length*UM_TO_MM,std_branch_length*UM_TO_MM);
    printf("[INFO] Total number of bifurcations = %u\n",angles.size());
    printf("[INFO] Bifurcation angle = %.2lf +/- %.2lf degrees\n",mean_biff_angle,std_biff_angle);

    write_geometric_info_to_file("outputs/geometric_info.txt",mean_segment_length,std_segment_length,\
                                mean_branch_length,std_branch_length,\
                                mean_biff_angle,std_biff_angle,\
                                segments.size(),branches.size(),angles.size());
}

void Graph::convert_to_initial_network_topology (const char filename[])
{
    std::vector< std::pair<uint32_t,uint32_t> > sorted_segments; 
    build_and_sort_segments(sorted_segments);
    
    std::vector<Point_Segment*> points;
    std::vector<Segment*> segments;
    build_shocker_structure_from_graph(sorted_segments,points,segments);
    
    write_shocker_initial_network_topology(filename,points,segments);
}

void Graph::get_segment_length (std::vector<double> &the_segments)
{
    uint32_t num_nodes = this->list_nodes.size();
    for (uint32_t i = 0; i < this->list_nodes.size(); i++)
    {
        uint32_t num_edges = this->list_nodes[i].list_edges.size();
        for (uint32_t j = 0; j < num_edges; j++)
        {
            double length = this->list_nodes[i].list_edges[j].length;
            the_segments.push_back(length);
        }
    }
}

void Graph::get_bifurcation_angles (std::vector<double> &the_angles)
{
    uint32_t num_nodes = this->list_nodes.size();
    for (uint32_t i = 0; i < this->list_nodes.size(); i++)
    {
        if (this->list_nodes[i].list_edges.size() == 2)
        {
            double u[3], v[3], angle;

            uint32_t p2_id = this->list_nodes[i].list_edges[0].dest_id;
            uint32_t p3_id = this->list_nodes[i].list_edges[1].dest_id;

            Node p1 = this->list_nodes[i];
            Node p2 = this->list_nodes[p2_id];
            Node p3 = this->list_nodes[p3_id];

            double a[3] = {p1.x,p1.y,p1.z};
            double b[3] = {p2.x,p2.y,p2.z};
            double c[3] = {p3.x,p3.y,p3.z};
            
            build_unitary_vector(u,a,b);
            build_unitary_vector(v,a,c);
            
            angle = calc_angle_between_vectors(u,v);

            the_angles.push_back(angle);
        }
    }
} 

void Graph::get_branches (std::vector<double> &the_branches)
{
    uint32_t num_branches = tag_branches();
    fill_branches(the_branches,num_branches);
    write_branches("outputs/branches.vtk");
}

uint32_t Graph::get_root_position_index (const double root_pos[])
{
    uint32_t root_id = 0;
    double closest_dist = __DBL_MAX__;
    for (uint32_t i = 0; i < this->list_nodes.size(); i++)
    {
        Node node = this->list_nodes[i];
        double dist = calc_norm(node.x,node.y,node.z,root_pos[0],root_pos[1],root_pos[2]);
        if (dist < closest_dist)
        {
            closest_dist = dist;
            root_id = i;
        }
    }
    return root_id;
}

double calculate_diameter (const double cv)
{
    static const double Gi = 7.9;
    static const double Cf = 3.4;
    static const double tauf = 0.1;
    double cv_m_per_s = cv / 1000.0;

    return (cv_m_per_s*cv_m_per_s*4.0*Cf*tauf) / (Gi) * 100.0;
}

double calculate_proportion (const double ref_diameter, const double ref_max_lat, const double aprox_max_lat)
{
    return ref_diameter * aprox_max_lat / ref_max_lat;
}

double calculate_propagation_velocity (const double diameter)
{
    static const double Gi = 7.9;
    static const double Cf = 3.4;
    static const double tauf = 0.1;

    return sqrt( (Gi*diameter)/(4.0*Cf*tauf) ) * 0.1;
}

double adjust_propagation_velocity (const double ref_cv, const double ref_max_lat, const double aprox_max_lat, double &new_diameter)
{
    double ref_diameter = calculate_diameter(ref_cv);
    new_diameter = calculate_proportion(ref_diameter,ref_max_lat,aprox_max_lat);
    return calculate_propagation_velocity(new_diameter)*1000.0;
}

void Graph::activate_edge (const uint32_t src_id, const uint32_t dest_id)
{
    // Search for the edge in the list of nodes
    uint32_t num_edges = this->list_nodes[src_id].list_edges.size();
    for (uint32_t i = 0; i < num_edges; i++)
    {
        if (this->list_nodes[src_id].list_edges[i].dest_id == dest_id)
        {
            this->list_nodes[src_id].list_edges[i].is_minimum_tree = true;
            return;
        }
    }
    printf("[-] ERROR! Edge (%u -> %u) was not found!\n",src_id,dest_id);
}

void Graph::fill_branches (std::vector<double> &the_branches, const uint32_t num_branches)
{
    for (uint32_t k = 0; k < num_branches; k++)
    {
        uint32_t cur_branch = k;
        double branch_length = 0.0;
        // Search for the edges that are part of the current branch
        for (uint32_t i = 0; i < this->total_nodes; i++)
        {
            uint32_t num_edges = this->list_nodes[i].list_edges.size();
            for (uint32_t j = 0; j < num_edges; j++)
            {
                uint32_t branch_id = this->list_nodes[i].list_edges[j].branch_id;
                double length = this->list_nodes[i].list_edges[j].length;
                if (branch_id == cur_branch)
                {
                    branch_length += length;
                }
            }
        }
        // Add the branch length to the output array
        the_branches.push_back(branch_length);
    }
}

void Graph::tag_minimum_network (std::vector<Point> pmjs)
{
    for (uint32_t i = 0; i < pmjs.size(); i++)
    {
        Point pmj = pmjs[i];
        uint32_t index = get_closest_terminal_point(pmj,false);
        Node tmp = this->list_nodes[index];
        while (this->parent[tmp.id] != -1)
        {
            uint32_t cur_id = tmp.id;
            uint32_t parent_id = this->parent[cur_id];
            
            activate_edge(parent_id,cur_id);
            tmp = this->list_nodes[parent_id];
        }
    }
}

uint32_t Graph::tag_branches ()
{
    uint32_t branch_counter = 1;
    for (uint32_t i = 0; i < this->total_nodes; i++)
    {
        uint32_t num_edges = this->list_nodes[i].list_edges.size();
        if (num_edges == 2)
        {
            // First branch starting from the bifurcation
            this->list_nodes[i].list_edges[0].branch_id = branch_counter;
            uint32_t next_id = this->list_nodes[i].list_edges[0].dest_id;
            Node *u = &this->list_nodes[next_id];
            while (u->list_edges.size() == 1)
            {
                next_id = u->list_edges[0].dest_id;
                u->list_edges[0].branch_id = branch_counter;
                u = &this->list_nodes[next_id];
            }
            branch_counter++;

            // Second branch starting from the bifurcation
            this->list_nodes[i].list_edges[1].branch_id = branch_counter;
            next_id = this->list_nodes[i].list_edges[1].dest_id;
            u = &this->list_nodes[next_id];
            while (u->list_edges.size() == 1)
            {
                next_id = u->list_edges[0].dest_id;
                u->list_edges[0].branch_id = branch_counter;
                u = &this->list_nodes[next_id];
            }
            branch_counter++;
        }
    }
    return branch_counter--;
}

uint32_t Graph::count_active_edges ()
{
    uint32_t num_active_edges = 0;
    for (uint32_t i = 0; i < this->list_nodes.size(); i++)
    {
        uint32_t num_edges = this->list_nodes[i].list_edges.size();
        for (uint32_t j = 0; j < num_edges; j++)
        {
            if (this->list_nodes[i].list_edges[j].is_minimum_tree)
                num_active_edges++;
        }
    }
    return num_active_edges;
}

void Graph::write_terminals (const char filename[])
{
    FILE *file = fopen(filename,"w+");
   
    fprintf(file,"# vtk DataFile Version 4.1\n");
    fprintf(file,"vtk output\n");
    fprintf(file,"ASCII\n");
    fprintf(file,"DATASET POLYDATA\n");
    
    uint32_t num_terminals = this->terminals_indexes.size();
    fprintf(file,"POINTS %u float\n",num_terminals);
    for (uint32_t i = 0; i < num_terminals; i++)
    {
        uint32_t id = this->terminals_indexes[i];
		fprintf(file,"%g %g %g\n",this->list_nodes[id].x,this->list_nodes[id].y,this->list_nodes[id].z);
    }
    
    fprintf(file,"VERTICES %u %u\n",num_terminals,num_terminals*2);
    for (uint32_t i = 0; i < num_terminals; i++)
        fprintf(file,"1 %u\n",i);
    
    fprintf(file,"POINT_DATA %u\n",num_terminals);
    fprintf(file,"FIELD FieldData 2\n");
    fprintf(file,"LAT 1 %u float\n",num_terminals);
    for (uint32_t i = 0; i < num_terminals; i++)
    {
        uint32_t id = this->terminals_indexes[i];
        if (this->lat[id] < this->min_LAT) this->min_LAT = this->lat[id];
		fprintf(file,"%g\n",this->lat[id]);
    }
    fprintf(file,"METADATA\n");
    fprintf(file,"INFORMATION 0\n\n");
    fprintf(file,"Active 1 %u float\n",num_terminals);
    for (uint32_t i = 0; i < num_terminals; i++)
    {
        fprintf(file,"1\n");
    }
    fprintf(file,"METADATA\n");
    fprintf(file,"INFORMATION 0\n");
    
    fclose(file);
}

void Graph::write_active_pmjs (const char filename[], const double perc)
{
    // The number of active PMJS is a percentage of the total number of terminals
    uint32_t num_terminals = this->terminals_indexes.size();
    const uint32_t num_pmjs = num_terminals*perc;           
    
    std::vector<bool> taken_points;
    taken_points.assign(num_terminals,false);

    FILE *file = fopen(filename,"w+");
    fprintf(file,"# vtk DataFile Version 4.1\n");
    fprintf(file,"vtk output\n");
    fprintf(file,"ASCII\n");
    fprintf(file,"DATASET POLYDATA\n");
    
    fprintf(file,"POINTS %u float\n",num_pmjs);

    uint32_t counter = 0;
    while (counter != num_pmjs)
    {
        uint32_t k = rand() % num_terminals;
        if (!taken_points[k])
        {
            uint32_t id = this->terminals_indexes[k];
		    fprintf(file,"%g %g %g\n",this->list_nodes[id].x,this->list_nodes[id].y,this->list_nodes[id].z);
            taken_points[k] = true;
            counter++;
        }    
    }
    
    fprintf(file,"VERTICES %u %u\n",num_pmjs,num_pmjs*2);
    for (uint32_t i = 0; i < num_pmjs; i++)
        fprintf(file,"1 %u\n",i);

    fclose(file);
}

void Graph::write_minimum_network (std::string pmj_filename, std::string filename)
{
    std::vector<Point> pmjs;
    read_points(pmj_filename,pmjs);

    depth_first_search_topology(this->root_index);

    tag_minimum_network(pmjs);
    uint32_t num_active_edges = count_active_edges();    

    // Write the minimum network
    // All the points are written, but only the minimum tree edges are considered
    FILE *file = fopen(filename.c_str(),"w+");
    fprintf(file,"# vtk DataFile Version 4.1\n");
    fprintf(file,"vtk output\n");
    fprintf(file,"ASCII\n");
    fprintf(file,"DATASET POLYDATA\n");
    fprintf(file,"POINTS %u float\n",this->list_nodes.size());
    for (uint32_t i = 0; i < this->list_nodes.size(); i++)
    {
        fprintf(file,"%g %g %g\n",this->list_nodes[i].x,this->list_nodes[i].y,this->list_nodes[i].z);
    }
    fprintf(file,"LINES %u %u\n",num_active_edges,num_active_edges*3);
    for (uint32_t i = 0; i < this->list_nodes.size(); i++)
    {
        uint32_t u = this->list_nodes[i].id;
        for (uint32_t j = 0; j < this->list_nodes[i].list_edges.size(); j++)
        {
            uint32_t v = this->list_nodes[i].list_edges[j].dest_id;
            if (this->list_nodes[i].list_edges[j].is_minimum_tree)
                fprintf(file,"2 %u %u\n",u,v);
        }
    }
    fclose(file);
}

void Graph::write_graph_info (const char filename[])
{
    depth_first_search_topology(this->root_index);

    FILE *file = fopen(filename,"w+");
    fprintf(file,"# vtk DataFile Version 4.1\n");
    fprintf(file,"vtk output\n");
    fprintf(file,"ASCII\n");
    fprintf(file,"DATASET POLYDATA\n");
    
    fprintf(file,"POINTS %u float\n",this->list_nodes.size());
    for (uint32_t i = 0; i < this->list_nodes.size(); i++)
    {
        fprintf(file,"%g %g %g\n",this->list_nodes[i].x,this->list_nodes[i].y,this->list_nodes[i].z);
    }
    fprintf(file,"LINES %u %u\n",this->total_edges,this->total_edges*3);
    for (uint32_t i = 0; i < this->list_nodes.size(); i++)
    {
        uint32_t u = this->list_nodes[i].id;
        for (uint32_t j = 0; j < this->list_nodes[i].list_edges.size(); j++)
        {
            uint32_t v = this->list_nodes[i].list_edges[j].dest_id;
            fprintf(file,"2 %u %u\n",u,v);
        }
    }
    fprintf(file,"POINT_DATA %u\n",this->list_nodes.size());
    fprintf(file,"FIELD FieldData 4\n");
    fprintf(file,"LAT 1 %u float\n",this->list_nodes.size());
    for (uint32_t i = 0; i < this->list_nodes.size(); i++)
        fprintf(file,"%g\n",this->lat[i]);
    fprintf(file,"METADATA\n");
    fprintf(file,"INFORMATION 0\n\n");
    fprintf(file,"Dist 1 %u float\n",this->list_nodes.size());
    for (uint32_t i = 0; i < this->list_nodes.size(); i++)
        fprintf(file,"%d\n",this->dist[i]);
    fprintf(file,"METADATA\n");
    fprintf(file,"INFORMATION 0\n\n");
    fprintf(file,"DFS_Counter 1 %u float\n",this->list_nodes.size());
    for (uint32_t i = 0; i < this->list_nodes.size(); i++)
        fprintf(file,"%d\n",dfs_counter[i]);
        fprintf(file,"Counter 1 %u float\n",this->list_nodes.size());
    for (uint32_t i = 0; i < this->list_nodes.size(); i++)
        fprintf(file,"%d\n",this->list_nodes[i].id);
    fprintf(file,"METADATA\n");
    fprintf(file,"INFORMATION 0\n\n");
    fclose(file);
}

void Graph::write_branches (const char filename[])
{
    FILE *file = fopen(filename,"w+");
    fprintf(file,"# vtk DataFile Version 4.1\n");
    fprintf(file,"vtk output\n");
    fprintf(file,"ASCII\n");
    fprintf(file,"DATASET POLYDATA\n");
    
    fprintf(file,"POINTS %u float\n",this->list_nodes.size());
    for (uint32_t i = 0; i < this->list_nodes.size(); i++)
    {
        fprintf(file,"%g %g %g\n",this->list_nodes[i].x,this->list_nodes[i].y,this->list_nodes[i].z);
    }
    fprintf(file,"LINES %u %u\n",this->total_edges,this->total_edges*3);
    for (uint32_t i = 0; i < this->list_nodes.size(); i++)
    {
        uint32_t u = this->list_nodes[i].id;
        for (uint32_t j = 0; j < this->list_nodes[i].list_edges.size(); j++)
        {
            uint32_t v = this->list_nodes[i].list_edges[j].dest_id;
            fprintf(file,"2 %u %u\n",u,v);
        }
    }
    fprintf(file,"CELL_DATA %u\n",this->total_edges);
    fprintf(file,"FIELD FieldData 1\n");
    fprintf(file,"Branch 1 %u float\n",this->total_edges);
    for (uint32_t i = 0; i < this->list_nodes.size(); i++)
    {
        uint32_t u = this->list_nodes[i].id;
        for (uint32_t j = 0; j < this->list_nodes[i].list_edges.size(); j++)
        {
            uint32_t branch_id = this->list_nodes[i].list_edges[j].branch_id;
            fprintf(file,"%u\n",branch_id);
        }
    }
    fclose(file);
}

void Graph::write_LAT (const char filename[])
{
    FILE *file = fopen(filename,"w+");
    fprintf(file,"# vtk DataFile Version 4.1\n");
    fprintf(file,"vtk output\n");
    fprintf(file,"ASCII\n");
    fprintf(file,"DATASET POLYDATA\n");
    
    fprintf(file,"POINTS %u float\n",this->list_nodes.size());
    for (uint32_t i = 0; i < this->list_nodes.size(); i++)
    {
        fprintf(file,"%g %g %g\n",this->list_nodes[i].x,this->list_nodes[i].y,this->list_nodes[i].z);
    }
    fprintf(file,"LINES %u %u\n",this->total_edges,this->total_edges*3);
    for (uint32_t i = 0; i < this->list_nodes.size(); i++)
    {
        uint32_t u = this->list_nodes[i].id;
        for (uint32_t j = 0; j < this->list_nodes[i].list_edges.size(); j++)
        {
            uint32_t v = this->list_nodes[i].list_edges[j].dest_id;
            fprintf(file,"2 %u %u\n",u,v);
        }
    }
    fprintf(file,"POINT_DATA %u\n",this->list_nodes.size());
    fprintf(file,"SCALARS LAT float\n");
    fprintf(file,"LOOKUP_TABLE default\n");
    for (uint32_t i = 0; i < this->list_nodes.size(); i++)
    {
        fprintf(file,"%g\n",this->lat[i]);
    }
    fclose(file);

    //for (uint32_t i = 0; i < this->terminals_indexes.size(); i++)
    //    printf("Active PMJ = %u || LAT = %g\n",this->terminals_indexes[i],this->lat[this->terminals_indexes[i]]);
}

void Graph::build_and_sort_segments (std::vector< std::pair<uint32_t,uint32_t> > &sorted_segments)
{
    for (uint32_t i = 0; i < this->total_nodes; i++)
    {
        uint32_t num_edges = this->list_nodes[i].list_edges.size();
        uint32_t u = this->list_nodes[i].id;
        for (uint32_t j = 0; j < num_edges; j++)
        {
            uint32_t v = this->list_nodes[i].list_edges[j].dest_id;
            std::pair<uint32_t,uint32_t> segment = std::make_pair(u,v);
            sorted_segments.push_back(segment);
        }
    }
    std::sort(sorted_segments.begin(),sorted_segments.end());
}

// TODO: Allow different diameter
void Graph::build_shocker_structure_from_graph(std::vector< std::pair<uint32_t,uint32_t> > sorted_segments,\
                                            std::vector<Point_Segment*> &points, std::vector<Segment*> &segments)
{
    const double diameter = 62.15; // {um}

    build_shocker_points(points);
    build_shocker_segments(sorted_segments,points,segments);

    update_diameter_and_LAT(segments,diameter);
    update_flags(points,segments);
}

void Graph::build_shocker_points (std::vector<Point_Segment*> &points)
{
    for (uint32_t i = 0; i < this->list_nodes.size(); i++)
    {
        double pos[3];
        pos[0] = this->list_nodes[i].x;
        pos[1] = this->list_nodes[i].y;
        pos[2] = this->list_nodes[i].z;
        Point_Segment *p = new Point_Segment(i,pos);
        points.push_back(p);
    }
}

void Graph::build_shocker_segments (std::vector< std::pair<uint32_t,uint32_t> > sorted_segments,\
                                std::vector<Point_Segment*> points, std::vector<Segment*> &segments)
{
    // Initialize the Shocker Segment array
    for (uint32_t i = 0; i < sorted_segments.size(); i++)
    {
        uint32_t src_id = sorted_segments[i].first;
        uint32_t dest_id = sorted_segments[i].second;
        Point_Segment *src = points[src_id];
        Point_Segment *dest = points[dest_id];

        Segment *s = new Segment(i,src,dest,NULL,NULL,NULL);
        segments.push_back(s);
    }

    // Set the 'left', 'right' and 'parent' pointers
    for (uint32_t i = 0; i < sorted_segments.size(); i++)
    {
        uint32_t cur_src = sorted_segments[i].first;
        uint32_t cur_dest = sorted_segments[i].second;

        // Find the 'offsprings' and 'parent' segment index
        std::vector<uint32_t> offsprings_segments;
        std::vector<uint32_t> parent_segments;
        for (uint32_t j = 0; j < sorted_segments.size(); j++)
        {
            if (sorted_segments[j].first == cur_dest)
                offsprings_segments.push_back(j);
            if (sorted_segments[j].second == cur_src)
                parent_segments.push_back(j);
        }

        // Set the 'offsprings' pointers by checking each case
        if (offsprings_segments.size() == 2)    // Case 1: Default the segment has two offsprings
        {
            uint32_t left_id = offsprings_segments[0];      // CONVENTION: The first offspring will be the 'left'
            uint32_t right_id = offsprings_segments[1];     // CONVENTION: The second offspring will be the 'right'

            // Update the 'left' and 'right' pointers for the current segment
            segments[i]->left = segments[left_id];
            segments[i]->right = segments[right_id];
        }
        else if (offsprings_segments.size() == 1)    // Case 2: The segment has one offspring
        {
            uint32_t left_id = offsprings_segments[0];      // CONVENTION: The first offspring will be the 'left'

            // Update the 'left' pointer for the current segment
            segments[i]->left = segments[left_id];
        }

        // Set the 'parent' pointer
        if (parent_segments.size() == 1)    
        {
            uint32_t parent_id = parent_segments[0];

            // Update the 'parent' pointer for the current segment
            segments[i]->parent = segments[parent_id];
        }
    }
}

void Graph::update_diameter_and_LAT(std::vector<Segment*> &segments, const double diameter)
{
    // Update the 'radius' from all segments
    for (uint32_t i = 0; i < segments.size(); i++)
    {
        Segment *cur_segment = segments[i];
        cur_segment->diameter = diameter;
    }
    // Update the 'LAT' from all segments
    for (uint32_t i = 0; i < segments.size(); i++)
    {
        Segment *cur_segment = segments[i];
        cur_segment->lat = cur_segment->calc_terminal_local_activation_time();
    }
}

void Graph::update_flags (std::vector<Point_Segment*> &points, std::vector<Segment*> &segments)
{
    for (uint32_t i = 0; i < segments.size(); i++)
    {
        Segment *cur_segment = segments[i];
        cur_segment->mintree = true;

        // All terminal distal point are set to 'active'
        if (cur_segment->is_terminal())
            cur_segment->dest->is_active = true;
    }
}

void Graph::write_shocker_initial_network_topology(const char filename[], std::vector<Point_Segment*> points, std::vector<Segment*> segments)
{
    uint32_t num_points = points.size();
    uint32_t num_segments = segments.size();

    FILE *file = fopen(filename,"w+");
    
    // Write points
    fprintf(file,"%u\n",num_points);
    for (uint32_t i = 0; i < num_points; i++)
        fprintf(file,"%u %g %g %g %d\n",points[i]->id,points[i]->pos[0],points[i]->pos[1],points[i]->pos[2],points[i]->is_active);

    // Write segments
    fprintf(file,"%u\n",num_segments);
    for (uint32_t i = 0; i < num_segments; i++)
    {
        uint32_t cur_seg_id = segments[i]->id;
        uint32_t src_id = segments[i]->src->id;
        uint32_t dest_id = segments[i]->dest->id;
        double lat = segments[i]->lat;
        double diameter = segments[i]->diameter;
        bool min_tree = segments[i]->mintree;
        fprintf(file,"%u %u %u %g %g %d\n",cur_seg_id,src_id,dest_id,lat,diameter,min_tree);
    }
    for (uint32_t i = 0; i < num_segments; i++)
    {
        int cur_id = segments[i]->id;
        int parent_id = (segments[i]->parent) ? segments[i]->parent->id : -1;
        int left_id = (segments[i]->left) ? segments[i]->left->id : -1;
        int right_id = (segments[i]->right) ? segments[i]->right->id : -1;
        fprintf(file,"%d %d\n",cur_id,parent_id);
        fprintf(file,"%d %d\n",cur_id,left_id);
        fprintf(file,"%d %d\n",cur_id,right_id);
    }
    fclose(file);
}

// This function will only work if a single geodesic pathway is given as input 
void Graph::write_polyline (const char filename[])
{
    FILE *file = fopen(filename,"w+");
    fprintf(file,"# vtk DataFile Version 4.1\n");
    fprintf(file,"vtk output\n");
    fprintf(file,"ASCII\n");
    fprintf(file,"DATASET POLYDATA\n");

    fprintf(file,"POINTS %u float\n",this->list_nodes.size());
    for (uint32_t i = 0; i < this->list_nodes.size(); i++)
    {
        fprintf(file,"%g %g %g\n",this->list_nodes[i].x,this->list_nodes[i].y,this->list_nodes[i].z);
    }
    fprintf(file,"LINES 1 %u\n",this->list_nodes.size()+1);
    fprintf(file,"%u ",this->list_nodes.size());
    for (uint32_t i = 0; i < this->list_nodes.size(); i++)
    {
        fprintf(file,"%u ",i);
    }
}