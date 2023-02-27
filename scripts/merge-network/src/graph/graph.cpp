#include "graph.h"

Graph::Graph ()
{

}

Graph::Graph (std::vector<Point> points, std::vector<Line> lines)
{
    uint32_t num_nodes = points.size();
    uint32_t num_edges = lines.size();

    this->list_nodes.assign(num_nodes,Node());
    for (uint32_t i = 0; i < num_nodes; i++)
    {
        uint32_t id = points[i].id;
        double *pos = points[i].pos;

        this->list_nodes[i].id = id;
        memcpy(this->list_nodes[i].pos,pos,sizeof(double)*3);
    }

    for (uint32_t i = 0; i < num_edges; i++)
    {
        uint32_t src = lines[i].src;
        uint32_t dest = lines[i].dest;
        double length = calc_norm(points[src].pos[0],points[src].pos[1],points[src].pos[2],\
                                points[dest].pos[0],points[dest].pos[1],points[dest].pos[2]);
        
        Edge e(dest,length);
        this->list_nodes[src].list_edges.push_back(e);
    }
    total_nodes = num_nodes;
    total_edges = num_edges;

    // Run a Dijkstra shortest path to fill the 'dist' and 'parent' arrays
    dijkstra(0);

    // Compute the LAT of all cells using the reference conduction velocity
    compute_activation_times(REF_CV);

    // Search for the terminal points and fill the 'terminal_indexes' array
    fill_terminal_indexes();

}

Graph::Graph (std::vector<Point> points, std::vector<Line> lines, std::vector<double> point_scalars)
{
    uint32_t num_nodes = points.size();
    uint32_t num_edges = lines.size();

    this->list_nodes.assign(num_nodes,Node());
    for (uint32_t i = 0; i < num_nodes; i++)
    {
        uint32_t id = points[i].id;
        double *pos = points[i].pos;
        double sigma = point_scalars[i];

        this->list_nodes[i].id = id;
        memcpy(this->list_nodes[i].pos,pos,sizeof(double)*3);
        this->list_nodes[i].sigma = sigma;
    }

    for (uint32_t i = 0; i < num_edges; i++)
    {
        uint32_t src = lines[i].src;
        uint32_t dest = lines[i].dest;
        double length = calc_norm(points[src].pos[0],points[src].pos[1],points[src].pos[2],\
                                points[dest].pos[0],points[dest].pos[1],points[dest].pos[2]);
        
        Edge e(dest,length);
        this->list_nodes[src].list_edges.push_back(e);
    }
    total_nodes = num_nodes;
    total_edges = num_edges;

    // Run a Dijkstra shortest path to fill the 'dist' and 'parent' arrays
    dijkstra(0);

    // Compute the LAT of all cells using the reference conduction velocity
    compute_activation_times(REF_CV);

    // Search for the terminal points and fill the 'terminal_indexes' array
    fill_terminal_indexes();

}

void Graph::depth_first_search (const uint32_t src_id, std::vector<double> &the_segments)
{
    std::vector<bool> dfs_visited;
    dfs_visited.assign(this->total_nodes,false);

    uint32_t total_segments = 0;
    double segment_size = 0.0;
    uint32_t flag = 0;

    dfs(this->list_nodes[src_id],dfs_visited,the_segments,segment_size,flag,total_segments);
}

void Graph::dfs (Node u, std::vector<bool> &dfs_visited, std::vector<double> &segments, double &segment_size, uint32_t &flag, uint32_t &total_segments)
{
    flag = 0;

    uint32_t u_index = u.id;
    dfs_visited[u_index] = true;

    uint32_t num_edges = u.list_edges.size();
    for (uint32_t j = 0; j < num_edges; j++)
    {
        uint32_t v_index = u.list_edges[j].dest_id;
        double length = u.list_edges[j].length;

        if (!dfs_visited[v_index])
        {
            flag++;
            segment_size += length;
            dfs(this->list_nodes[v_index],dfs_visited,segments,segment_size,flag,total_segments);
        }
    }
    // Terminal edge
    if (u.list_edges.size() == 0)
    {
        total_segments++;
        segments.push_back(segment_size);
        segment_size = 0.0;
    }
    if (flag == 2)
    {
        total_segments++;
    }
}

void Graph::print ()
{
    for (uint32_t i = 0; i < this->total_nodes; i++)
    {
        this->list_nodes[i].print();
    }
}

uint32_t Graph::get_closest_point (Point p)
{
    uint32_t closest_index = 0;
    double closest_dist = __DBL_MAX__;
    for (uint32_t i = 0; i < this->list_nodes.size(); i++)
    {
        Node node = this->list_nodes[i];
        double dist = calc_norm(node.pos[0],node.pos[1],node.pos[2],p.pos[0],p.pos[1],p.pos[2]);
        if (dist < closest_dist && this->list_nodes[i].list_edges.size() > 0)
        {
            closest_dist = dist;
            closest_index = i;
        }
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

void Graph::fill_terminal_indexes ()
{
    uint32_t num_nodes = this->list_nodes.size();
    for (uint32_t i = 0; i < num_nodes; i++)
    {
        if (this->list_nodes[i].list_edges.size() == 0)
            this->terminals_indexes.push_back(i);
    }
}

void Graph::compute_geometrics ()
{
    std::vector<double> segments;
    get_segment_length(segments);
    write_data_to_file("outputs/all_segment_length.dat",segments);

    std::vector<double> angles;
    get_bifurcation_angles(angles);
    write_data_to_file("outputs/all_bifurcation_angles.dat",angles);

    double mean_segment_length, std_segment_length;
    compute_mean_std(segments,mean_segment_length,std_segment_length);

    double mean_biff_angle, std_biff_angle;
    compute_mean_std(angles,mean_biff_angle,std_biff_angle);

    printf("[INFO] Total number of segment = %u\n",segments.size());
    printf("[INFO] Segment length = %.2lf +/- %.2lf mm\n",mean_segment_length*UM_TO_MM,std_segment_length*UM_TO_MM);
    printf("[INFO] Total number of bifurcations = %u\n",angles.size());
    printf("[INFO] Bifurcation angle = %.2lf +/- %.2lf degrees\n",mean_biff_angle,std_biff_angle);
}

void Graph::compute_activation_times (const double cv)
{
    this->lat.clear();
    for (uint32_t i = 0; i < this->total_nodes; i++)
    {
        if (this->dist[i] == __DBL_MAX__)
            this->lat.push_back(0);
        else
            this->lat.push_back(this->dist[i]/cv);
    }
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

            double a[3] = {p1.pos[0],p1.pos[1],p1.pos[2]};
            double b[3] = {p2.pos[0],p2.pos[1],p2.pos[2]};
            double c[3] = {p3.pos[0],p3.pos[1],p3.pos[2]};
            
            build_unitary_vector(u,a,b);
            build_unitary_vector(v,a,c);
            
            angle = calc_angle_between_vectors(u,v);

            the_angles.push_back(angle);
        }
    }
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
		fprintf(file,"%g %g %g\n",this->list_nodes[id].pos[0],this->list_nodes[id].pos[1],this->list_nodes[id].pos[2]);
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

void Graph::write_network (const char filename[])
{
    FILE *file = fopen(filename,"w+");
    fprintf(file,"# vtk DataFile Version 4.1\n");
    fprintf(file,"vtk output\n");
    fprintf(file,"ASCII\n");
    fprintf(file,"DATASET POLYDATA\n");
    
    fprintf(file,"POINTS %u float\n",this->list_nodes.size());
    for (uint32_t i = 0; i < this->list_nodes.size(); i++)
        fprintf(file,"%g %g %g\n",this->list_nodes[i].pos[0],this->list_nodes[i].pos[1],this->list_nodes[i].pos[2]);
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
        fprintf(file,"%g %g %g\n",this->list_nodes[i].pos[0],this->list_nodes[i].pos[1],this->list_nodes[i].pos[2]);
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
        fprintf(file,"%g\n",this->lat[i]);
    fclose(file);
}

void Graph::write_MonoAlg3D (const char filename[])
{
    FILE *file = fopen(filename,"w+");
    fprintf(file,"# vtk DataFile Version 4.1\n");
    fprintf(file,"vtk output\n");
    fprintf(file,"ASCII\n");
    fprintf(file,"DATASET POLYDATA\n");
    
    fprintf(file,"POINTS %u float\n",this->list_nodes.size());
    for (uint32_t i = 0; i < this->list_nodes.size(); i++)
        fprintf(file,"%g %g %g\n",this->list_nodes[i].pos[0],this->list_nodes[i].pos[1],this->list_nodes[i].pos[2]);
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
    fprintf(file,"SCALARS sigma float\n");
    fprintf(file,"LOOKUP_TABLE default\n");
    for (uint32_t i = 0; i < this->list_nodes.size(); i++)
        fprintf(file,"%g\n",this->list_nodes[i].sigma);
    fclose(file);
}

Graph* merge_networks (Graph *his, Graph *lv, Graph *rv, std::vector<Point> the_roots)
{
    uint32_t num_nodes_his = his->list_nodes.size();        // Eliminate the root points of the LV and RV in the His-bundle
    uint32_t num_nodes_lv = lv->list_nodes.size();
    uint32_t num_nodes_rv = rv->list_nodes.size();
    uint32_t num_nodes_merged = num_nodes_his + num_nodes_lv + num_nodes_rv;
    uint32_t num_edges_merged = his->total_edges + lv->total_edges + rv->total_edges; 

    Graph *merged_network = new Graph();
    merged_network->list_nodes.assign(num_nodes_merged,Node());
    merged_network->total_nodes = num_nodes_merged;
    
    // His-Bundle 
    uint32_t offset = 0;
    for (uint32_t i = 0; i < num_nodes_his; i++)
    {
        uint32_t id = i + offset;
        double *pos = his->list_nodes[i].pos;
        double sigma = his->list_nodes[i].sigma;
        merged_network->list_nodes[i].setNode(id,pos,sigma);

        uint32_t num_edges = his->list_nodes[i].list_edges.size();
        uint32_t src = i;
        for (uint32_t j = 0; j < num_edges; j++)
        {
            uint32_t dest = his->list_nodes[i].list_edges[j].dest_id + offset;
            double length = calc_norm(his->list_nodes[src].pos[0],his->list_nodes[src].pos[1],his->list_nodes[src].pos[2],\
                                his->list_nodes[dest].pos[0],his->list_nodes[dest].pos[1],his->list_nodes[dest].pos[2]);
            Edge e(dest,length);
            merged_network->list_nodes[id].list_edges.push_back(e);
        }
    }

    // LV
    offset = num_nodes_his;
    for (uint32_t i = 0; i < num_nodes_lv; i++)
    {
        uint32_t id = i + offset;
        double *pos = lv->list_nodes[i].pos;
        double sigma = lv->list_nodes[i].sigma;
        merged_network->list_nodes[id].setNode(id,pos,sigma);

        uint32_t num_edges = lv->list_nodes[i].list_edges.size();
        uint32_t src_lv = i;
        for (uint32_t j = 0; j < num_edges; j++)
        {
            uint32_t dest_lv = lv->list_nodes[i].list_edges[j].dest_id; 
            double length = calc_norm(lv->list_nodes[src_lv].pos[0],lv->list_nodes[src_lv].pos[1],lv->list_nodes[src_lv].pos[2],\
                                lv->list_nodes[dest_lv].pos[0],lv->list_nodes[dest_lv].pos[1],lv->list_nodes[dest_lv].pos[2]);
            uint32_t dest = dest_lv + offset;
            Edge e(dest,length);
            merged_network->list_nodes[id].list_edges.push_back(e);
        }
    }
    
    // RV
    offset = num_nodes_his + num_nodes_lv;
    for (uint32_t i = 0; i < num_nodes_rv; i++)
    {
        uint32_t id = i + offset;
        double *pos = rv->list_nodes[i].pos;
        double sigma = rv->list_nodes[i].sigma;
        merged_network->list_nodes[id].setNode(id,pos,sigma);

        uint32_t num_edges = rv->list_nodes[i].list_edges.size();
        uint32_t src_rv = i;
        for (uint32_t j = 0; j < num_edges; j++)
        {
            uint32_t dest_rv = rv->list_nodes[i].list_edges[j].dest_id;
            double length = calc_norm(rv->list_nodes[src_rv].pos[0],rv->list_nodes[src_rv].pos[1],rv->list_nodes[src_rv].pos[2],\
                                rv->list_nodes[dest_rv].pos[0],rv->list_nodes[dest_rv].pos[1],rv->list_nodes[dest_rv].pos[2]);
            uint32_t dest = dest_rv + offset;
            Edge e(dest,length);
            merged_network->list_nodes[id].list_edges.push_back(e);
        }
    }

    // Link the LV and RV Purkinje networks
    uint32_t lv_root_index = merged_network->get_closest_point(the_roots[0]);
    uint32_t rv_root_index = merged_network->get_closest_point(the_roots[1]);

    double length;
    length = calc_norm(merged_network->list_nodes[1].pos[0],merged_network->list_nodes[1].pos[1],merged_network->list_nodes[1].pos[2],\
                    merged_network->list_nodes[lv_root_index].pos[0],merged_network->list_nodes[lv_root_index].pos[1],merged_network->list_nodes[lv_root_index].pos[2]);
    Edge his_lv(lv_root_index,length);
    merged_network->list_nodes[1].list_edges.push_back(his_lv);

    length = calc_norm(merged_network->list_nodes[1].pos[0],merged_network->list_nodes[1].pos[1],merged_network->list_nodes[1].pos[2],\
                    merged_network->list_nodes[rv_root_index].pos[0],merged_network->list_nodes[rv_root_index].pos[1],merged_network->list_nodes[rv_root_index].pos[2]);
    Edge his_rv(rv_root_index,length);
    merged_network->list_nodes[1].list_edges.push_back(his_rv);

    merged_network->total_nodes = num_nodes_merged;
    merged_network->total_edges = num_edges_merged + 2;

    // Run a Dijkstra shortest path to fill the 'dist' and 'parent' arrays
    merged_network->dijkstra(0);

    // Compute the LAT of all cells using the reference conduction velocity
    merged_network->compute_activation_times(merged_network->REF_CV);

    // Search for the terminal points and fill the 'terminal_indexes' array
    merged_network->fill_terminal_indexes();

    return merged_network;

}