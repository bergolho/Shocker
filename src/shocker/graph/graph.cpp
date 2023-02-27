#include "graph.h"

GraphNode::GraphNode (const u_int32_t id, const double pos[])
{
    this->id = id;
    memcpy(this->pos,pos,sizeof(double)*3);
    this->num_edges = 0;
    this->next = nullptr;
    this->list_edges = nullptr;
}

GraphNode::~GraphNode ()
{
    this->list_edges = nullptr;
    this->next = nullptr;
}

void GraphNode::print ()
{
    printf("Node %u -> ( %g %g %g ) {%u} ",this->id,this->pos[0],this->pos[1],this->pos[2],this->num_edges);
}

GraphEdge::GraphEdge (const u_int32_t id, const double length, GraphNode *dest)
{
    this->id = id;
    this->length = length;
    this->dest = dest;
    this->next = nullptr;
}

GraphEdge::~GraphEdge ()
{
    this->dest = nullptr;
    this->next = nullptr;
}

void GraphEdge::print ()
{
    printf("Edge %u -> length = %g",this->id,this->length);
}

Graph::Graph ()
{
    this->last_node = nullptr;
    this->list_nodes = nullptr;
    this->total_nodes = 0;
    this->total_edges = 0;
}

Graph::~Graph ()
{
    if (this->list_nodes)
        free_list_nodes(this->list_nodes);
}

void Graph::free_list_nodes (GraphNode *head)
{
    GraphNode *n1 = head;
    GraphNode *n2 = head->next;
    while (n1 != nullptr)
    {
        if (n1->list_edges)
            free_list_edges(n1);
        free(n1);
        n1 = n2;
        if (n2 != nullptr)
            n2 = n2->next;
    }
}

void Graph::free_list_edges (GraphNode *node)
{
    GraphEdge *e1 = node->list_edges;
    GraphEdge *e2 = node->list_edges->next;
    while (e1 != nullptr)
    {
        free(e1);
        e1 = e2;
        if (e2 != nullptr)
            e2 = e2->next;
    }
    node->list_edges = nullptr;
}

void Graph::insert_node (const double pos[])
{
    GraphNode *tmp = this->list_nodes;
    GraphNode *node = new GraphNode(this->total_nodes++,pos);
    if (!tmp)
    {
        this->list_nodes = node;
        this->last_node = node;
    }
    else
    {
        this->last_node->next = node;
        this->last_node = this->last_node->next;
    }
}

void Graph::insert_edge (const uint32_t src_id, const uint32_t dst_id)
{
    // Check if the edge is invalid
    if (src_id == dst_id) return;

    GraphNode *n1 = search_node(src_id);
    GraphNode *n2 = search_node(dst_id);

    double length = euclidean_norm(n1->pos[0],n1->pos[1],n1->pos[2],n2->pos[0],n2->pos[1],n2->pos[2]);
    GraphEdge *edge = new GraphEdge(dst_id,length,n2);
    if (!n1->list_edges)
    {
        n1->list_edges = edge;
    }
    else
    {
        GraphEdge *e = n1->list_edges;
        while (e->next != nullptr)
            e = e->next;
        e->next = edge;
    }
    n1->num_edges++;
    this->total_edges++;
}

GraphNode* Graph::search_node (const uint32_t id)
{
    GraphNode *tmp = this->list_nodes;
    while (tmp != nullptr)
    {
        if (tmp->id == id) return tmp;
        tmp = tmp->next;
    }
    fprintf(stderr,"[graph] ERROR! Node %d was not found!\n",id);
    return nullptr;
}

Graph* convert_shocker_to_graph (std::vector<Node*> node_list, std::vector<Segment*> segment_list)
{
    Graph *result = new Graph();
    for (uint32_t i = 0; i < node_list.size(); i++)
        result->insert_node(node_list[i]->pos);
    for (uint32_t i = 0; i < segment_list.size(); i++)
    {
        result->insert_edge(segment_list[i]->src->id,segment_list[i]->dest->id);
    }
    return result;
}

Graph* build_mesh_purkinje (Graph *network, const double dx)
{
    Graph *result = new Graph();
    result->dx = dx;

    uint32_t n = network->total_nodes;
    uint32_t *map_network_to_mesh = (uint32_t*)calloc(n,sizeof(uint32_t));
    bool *dfs_visited = (bool*)malloc(sizeof(bool)*n);
    for (uint32_t i = 0; i < n; i++) dfs_visited[i] = false;

    GraphNode *tmp = network->list_nodes;
    result->insert_node(tmp->pos);
    result->depth_first_search(tmp,0,map_network_to_mesh,dfs_visited);

    return result;
}

void Graph::depth_first_search (GraphNode *u, const uint32_t level, uint32_t *map_network_to_mesh, bool *dfs_visited)
{
    dfs_visited[u->id] = true;
    GraphEdge *v = u->list_edges;
    while (v != nullptr)
    {
        if (dfs_visited[v->id] == false)
        {
            this->grow_segment(u,v,map_network_to_mesh);
            depth_first_search(v->dest,level+1,map_network_to_mesh,dfs_visited);
        }
        v = v->next;
    }
}

void Graph::grow_segment (GraphNode *u, GraphEdge *v, uint32_t *map_network_to_mesh) 
{
    double d_ori[3], d[3];
    double segment_length = v->length;
    double remainder_points = fmod(segment_length,this->dx);
    uint32_t n_points = segment_length / dx;

    // Capture the index of the growing node on the mesh
    uint32_t id_source = map_network_to_mesh[u->id];

    // Calculate a unitary direction vector of the segment
    calc_unitary_vector(d_ori,u,v->dest);

    // Copy the position of the source node
    memcpy(d,u->pos,sizeof(double)*3);

    // Grow the number of points of size 'h' until reaches the size of the segment
    for (uint32_t k = 1; k <= n_points; k++) {

        double pos[3];
        pos[0] = d[0] + d_ori[0]*dx*k;
        pos[1] = d[1] + d_ori[1]*dx*k;
        pos[2] = d[2] + d_ori[2]*dx*k;

        this->insert_node(pos);
        this->insert_edge(id_source,this->total_nodes-1);

        id_source = this->total_nodes-1;
    }
    // Grow any remainder point with a size 'dx'
    if (remainder_points > 0.0) {

        double pos[3];
        pos[0] = d[0] + d_ori[0]*dx*n_points + d_ori[0]*remainder_points;
        pos[1] = d[1] + d_ori[1]*dx*n_points + d_ori[1]*remainder_points;
        pos[2] = d[2] + d_ori[2]*dx*n_points + d_ori[2]*remainder_points;

        this->insert_node(pos);
        this->insert_edge(id_source,this->total_nodes-1);

        id_source = this->total_nodes-1;
    }

    // Save the last inserted node index, in case this node generates offsprings
    map_network_to_mesh[v->id] = id_source;
}

void calc_unitary_vector (double d_ori[], GraphNode *u, GraphNode *v)
{
    for (uint32_t i = 0; i < 3; i++) d_ori[i] = v->pos[i] - u->pos[i];
    double norm = sqrt(d_ori[0]*d_ori[0] + d_ori[1]*d_ori[1] + d_ori[2]*d_ori[2]);
    for (uint32_t i = 0; i < 3; i++) d_ori[i] /= norm;
}

void Graph::print ()
{
    GraphNode *n = this->list_nodes;
    while (n != nullptr)
    {
        GraphEdge *e = n->list_edges;
        fprintf(stdout,"|| %u ( %g %g %g ) ||",n->id,n->pos[0],n->pos[1],n->pos[2]);
        while (e != nullptr)
        {
            fprintf(stdout," --> || %u %g ( %g %g %g ) ||",e->id,e->length,e->dest->pos[0],e->dest->pos[1],e->dest->pos[2]);
            e = e->next;
        }
        fprintf(stdout,"\n");
        n = n->next;
    }
    printf("Total number of nodes = %u\n",this->total_nodes);
    printf("Total number of edges = %u\n",this->total_edges);
}

void Graph::write (std::string output_dir)
{
    std::string filename = output_dir + "/purkinje_mesh.vtk";
    FILE *file = fopen(filename.c_str(),"w+");
    fprintf(file,"# vtk DataFile Version 3.0\n");
    fprintf(file,"Purkinje\nASCII\nDATASET POLYDATA\n");
    fprintf(file,"POINTS %d float\n",this->total_nodes);

    GraphNode *n = this->list_nodes;
    while (n != nullptr)
    {
        fprintf(file,"%g %g %g\n",n->pos[0],n->pos[1],n->pos[2]);
        n = n->next;
    }
    n = this->list_nodes;
    fprintf(file,"LINES %u %u\n",this->total_edges,this->total_edges*3);
    while (n != nullptr)
    {
        GraphEdge *e = n->list_edges;
        while (e != nullptr)
        {
            fprintf(file,"2 %u %u\n",n->id,e->id);
            e = e->next;
        }
        n = n->next;
    }
    fclose(file);
}

void Graph::build_and_sort_segments (std::vector< std::pair<uint32_t,uint32_t> > &sorted_segments)
{
    GraphNode *n = this->list_nodes;
    while (n != nullptr)
    {
        uint32_t num_edges = n->num_edges;
        uint32_t u = n->id;
        GraphEdge *e = n->list_edges;
        while (e != nullptr)
        {
            uint32_t v = e->id;
            std::pair<uint32_t,uint32_t> segment = std::make_pair(u,v);
            sorted_segments.push_back(segment);
            e = e->next;
        }
        n = n->next;
    }
    std::sort(sorted_segments.begin(),sorted_segments.end());
}

void Graph::build_shocker_structure_from_graph(std::vector< std::pair<uint32_t,uint32_t> > sorted_segments,\
                                            const double diameter, const double lat_offset,\
                                            std::vector<Node*> &p_nodes, std::vector<Segment*> &p_segments)
{
    build_shocker_nodes(p_nodes);
    build_shocker_segments(sorted_segments,diameter,p_nodes,p_segments);

    update_diameter_and_LAT(p_segments,diameter,lat_offset);
    update_flags(p_nodes,p_segments);
}

void Graph::build_shocker_nodes (std::vector<Node*> &nodes)
{
    GraphNode *n = this->list_nodes;
    while (n != nullptr)
    {
        Node *node = new Node(n->id,n->pos);
        nodes.push_back(node);
        n = n->next;
    }
}

void Graph::build_shocker_segments (std::vector< std::pair<uint32_t,uint32_t> > sorted_segments, const double diameter,\
                                std::vector<Node*> nodes, std::vector<Segment*> &segments)
{
    // Initialize the Shocker Segment array
    for (uint32_t i = 0; i < sorted_segments.size(); i++)
    {
        uint32_t src_id = sorted_segments[i].first;
        uint32_t dest_id = sorted_segments[i].second;
        Node *src = nodes[src_id];
        Node *dest = nodes[dest_id];

        Segment *s = new Segment(i,diameter,0.0,true,src,dest,nullptr,nullptr,nullptr);
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

void Graph::update_diameter_and_LAT(std::vector<Segment*> &segments, const double diameter, const double lat_offset)
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
        cur_segment->lat = cur_segment->calc_terminal_local_activation_time() + lat_offset;
    }
}

void Graph::update_flags (std::vector<Node*> &nodes, std::vector<Segment*> &segments)
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