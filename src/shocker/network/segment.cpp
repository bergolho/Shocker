#include "segment.h"

bool first_call = true;
struct stop_watch calc_local_activation_terminal_time;
uint32_t total_time_local_activation_terminal = 0;

Segment::Segment ()
{
    this->id = 0;
    this->src = nullptr;
    this->dest = nullptr;
    this->parent = nullptr;
    this->left = nullptr;
    this->right = nullptr;
    this->diameter = 0;
    this->length = 0;
    this->lat = 0;
    this->mintree = false;
}

Segment::Segment (const uint32_t id, const double diameter, const double lat, const bool mintree,\
                Node *src, Node *dest,\
                Segment *left, Segment *right, Segment *parent)
{
    this->id = id;
    this->src = src;
    this->length = euclidean_norm(src->pos[0],src->pos[1],src->pos[2],dest->pos[0],dest->pos[1],dest->pos[2]);
    this->dest = dest;
    this->left = left;
    this->right = right;
    this->parent = parent;
    this->diameter = diameter;
    this->lat = lat;
    this->mintree = mintree;
    calc_middle_point(this->middle_pos);
}

Segment::Segment (Segment *in)
{
    this->id = in->id;
    this->src = in->src;
    this->dest = in->dest;
    this->length = in->length;
    this->diameter = in->diameter;
    this->lat = in->lat;
    this->mintree = in->mintree;
    this->left = in->left;
    this->right = in->right;
    this->parent = in->parent;
}

Segment::~Segment ()
{
    this->parent = nullptr;
    this->left = nullptr;
    this->right = nullptr;
    this->src = nullptr;
    this->dest = nullptr;
}

void Segment::calc_middle_point (double pos[])
{
    Node *src = this->src;
    Node *dest = this->dest;

    pos[0] = (src->pos[0] + dest->pos[0]) / 2.0;
    pos[1] = (src->pos[1] + dest->pos[1]) / 2.0;
    pos[2] = (src->pos[2] + dest->pos[2]) / 2.0;
}

void Segment::calc_unitary_vector (double u[])
{
    if (this->length < 1.0e-08) this->length = 1.0e-08;

    u[0] = (this->dest->pos[0] - this->src->pos[0]) / this->length;
    u[1] = (this->dest->pos[1] - this->src->pos[1]) / this->length;
    u[2] = (this->dest->pos[2] - this->src->pos[2]) / this->length;
}

void Segment::tag_pathway ()
{
    Segment *tmp = this;
    while (tmp != nullptr)
    {
        tmp->mintree = true;
        tmp = tmp->parent;
    }
}

double Segment::calc_radius ()
{
    return this->diameter/2.0;
}

double Segment::calc_length ()
{
    return euclidean_norm(src->pos[0],src->pos[1],src->pos[2],dest->pos[0],dest->pos[1],dest->pos[2]);
}

double Segment::calc_propagation_velocity ()
{
    const double G = 7.9;
    const double Cf = 3.4;
    const double tauf = 0.1;

    double d = this->diameter; 

    // Output in {m/s}
    return pow( (G*d)/(4.0*Cf*tauf) , 0.5 ) * 0.1;
}

// Input: {m/s}
// Output: {um} 
double Segment::calc_diameter (const double cv)
{
    const double G = 7.9;
    const double Cf = 3.4;
    const double tauf = 0.1;
    double v = cv*10.0;

    // Output in {um}
    return (4.0*v*v*Cf*tauf/G);
}

// Input: {m/s}
// Output: {um}
// Calculate the corresponding 'sigma' for the diameter for a MonoAlg3D simulation
// The values were fitted linearly using the Trovato_2019 model and a 10cm cable
double Segment::calc_conductivity ()
{
    const double a = 2.97e-05;
    const double b = -1.85e-04;
    return a*this->diameter+b;
}

double Segment::calc_local_activation_time ()
{
    double length = this->length;                           // {um}
    double cv = calc_propagation_velocity()*M_S_TO_UM_MS;   // {um/ms}
    
    // Output in {ms}
    return length / cv;
}

double Segment::calc_pathway_length ()
{
    double result = 0.0;
    Segment *tmp = this;
    while (tmp != nullptr)
    {
        result += tmp->length;
        tmp = tmp->parent;
    }
    return result;
}

double Segment::calc_pathway_length (std::vector<double> &out)
{
    double result = 0.0;
    Segment *tmp = this;
    while (tmp != nullptr)
    {
        result += tmp->length;
        out.push_back(tmp->length);
        tmp = tmp->parent;
    }
    return result;
}

double Segment::calc_terminal_local_activation_time ()
{
    start_stop_watch(&calc_local_activation_terminal_time);

    double result = 0.0;
    Segment *tmp = this;
    while (tmp != nullptr)
    {
        double lat = tmp->calc_local_activation_time();
        result += lat;
        tmp = tmp->parent;
    }

    total_time_local_activation_terminal += stop_stop_watch(&calc_local_activation_terminal_time);

    // Return the LAT in {ms}
    return result; 
}

double Segment::calc_dproj (const double pos[])
{
    Node *prox = this->src;
    Node *dist = this->dest;
    double length = this->length;

    double dot_product = (prox->pos[0] - dist->pos[0])*(pos[0] - dist->pos[0]) +\
                         (prox->pos[1] - dist->pos[1])*(pos[1] - dist->pos[1]) +\
                         (prox->pos[2] - dist->pos[2])*(pos[2] - dist->pos[2]);

    return dot_product * pow(length,-2.0);
}

double Segment::calc_dortho (const double pos[])
{
    Node *prox = this->src;
    Node *dist = this->dest;
    double length = this->length;

    double length_term = euclidean_norm(pos[0],pos[1],pos[2],\
                                dist->pos[0],dist->pos[1],dist->pos[2]);

    double dot_product = (prox->pos[0] - dist->pos[0])*(dist->pos[0] - pos[0]) +\
                         (prox->pos[1] - dist->pos[1])*(dist->pos[1] - pos[1]) +\
                         (prox->pos[2] - dist->pos[2])*(dist->pos[2] - pos[2]);

    return sqrt( powf(length_term * length, 2.0) - powf(dot_product,2.0) ) * powf(length,-1.0);
}

double Segment::calc_dend (const double pos[])
{
    Node *prox = this->src;
    Node *dist = this->dest;
    double length = this->length;

    double d_dist = euclidean_norm(pos[0],pos[1],pos[2],\
                                dist->pos[0],dist->pos[1],dist->pos[2]);

    double d_prox = euclidean_norm(pos[0],pos[1],pos[2],\
                                prox->pos[0],prox->pos[1],prox->pos[2]);

    return std::min(d_dist,d_prox);
}

uint32_t Segment::calc_level ()
{
    uint32_t result = 0;
    Segment *tmp = this;
    while (tmp != nullptr)
    {
        result++;
        tmp = tmp->parent;
    }
    return result;
}

bool Segment::is_terminal ()
{
    return (this->left == nullptr && this->right == nullptr) ? true : false;
}

bool Segment::is_bifurcation ()
{
    return (this->left != nullptr && this->right != nullptr) ? true : false;
}

bool Segment::is_segment ()
{
    return (this->right != nullptr && this->left == nullptr) ? true : false;
}

bool Segment::is_mintree ()
{
    return this->mintree;
}

bool Segment::is_inside_region (const double center[], const double radius)
{
    double dist = euclidean_norm(dest->pos[0],dest->pos[1],dest->pos[2],center[0],center[1],center[2]); 
    return (dist < radius) ? true : false;
}

void Segment::print ()
{
    printf("Segment %u (%u,%u) (%g %g %g) - (%g %g %g) -- DIAMETER = %g μm -- LAT = %g ms\n",this->id,src->id,dest->id,\
                                                                        src->pos[0],src->pos[1],src->pos[2],dest->pos[0],dest->pos[1],dest->pos[2],\
                                                                        this->diameter,this->lat);

    if (this->parent == nullptr)
        printf("\tPARENT = NIL\n");
    else
        printf("\tPARENT = %u\n",this->parent->id);
    if (this->left == nullptr)
        printf("\tLEFT = NIL\n");
    else
        printf("\tLEFT = %u\n",this->left->id);
    if (this->right == nullptr)
        printf("\tRIGHT = NIL\n");
    else
        printf("\tRIGHT = %u\n",this->right->id);
}

void eliminate_segment_from_list (std::vector<Segment*> &s_list, Segment *s)
{
    s_list.erase(s_list.begin() + s->id);
    order_segment_list(s_list);
    delete s;
    s = nullptr;
}

void order_segment_list (std::vector<Segment*> &s_list)
{
    uint32_t counter = 0;
    for (uint32_t i = 0; i < s_list.size(); i++)
    {
        s_list[i]->id = counter;
        counter++;
    }
}

void get_segment_length (std::vector<Segment*> s_list, std::vector<double> &segments)
{
    uint32_t ns = s_list.size();
    for (uint32_t i = 0; i < ns; i++)
    {
        Segment *s = s_list[i];
        segments.push_back(s->length*UM_TO_MM);    // Out: {mm}
    }
}

void get_bifurcation_angles(std::vector<Segment*> s_list, std::vector<double> &angles)
{
    uint32_t ns = s_list.size();
    for (uint32_t i = 0; i < ns; i++)
    {
        Segment *s = s_list[i];
        if (s->is_bifurcation())
        {
            double u[3], v[3], angle;
            s->left->calc_unitary_vector(u);
            s->right->calc_unitary_vector(v);
            angle = calc_angle_between_vectors(u,v);
            angles.push_back(angle);               // Out: {degrees}
        }
    }
}

Segment* get_closest_segment_by_distance (std::vector<Segment*> s_list, const double pos[])
{
    uint32_t ns = s_list.size();
    double middle_pos[3];
    double min_dist = __DBL_MAX__;
    Segment *closest_segment = nullptr;
    for (uint32_t i = 0; i < ns; i++)
    {
        Segment *s = s_list[i];
        s->calc_middle_point(middle_pos);
        double dist = euclidean_norm(middle_pos[0],middle_pos[1],middle_pos[2],\
                                    pos[0],pos[1],pos[2]);
        if (dist < min_dist)
        {
            min_dist = dist;
            closest_segment = s;
        }
    }
    return closest_segment;
}

std::vector<Segment*> get_closest_segment_list_by_distance (std::vector<Segment*> s_list, const double pos[], const uint32_t num_size)
{
    std::vector<Segment*> result;
    uint32_t ns = s_list.size();
    double middle_pos[3];
    double min_dist = __DBL_MAX__;
    std::vector<std::pair<double,uint32_t>> dist_table;
    for (uint32_t i = 0; i < ns; i++)
    {
        Segment *s = s_list[i];
        s->calc_middle_point(middle_pos);
        double dist = euclidean_norm(middle_pos[0],middle_pos[1],middle_pos[2],pos[0],pos[1],pos[2]);
        dist_table.push_back(std::make_pair(dist,i));
    }
    std::sort(dist_table.begin(),dist_table.end());
    //for (uint32_t i = 0; i < dist_table.size(); i++) // All segments
    // Pick up only the first "N_a" segments
    for (uint32_t i = 0; i < num_size; i++)
    {
        uint32_t id = dist_table[i].second;
        Segment *s = s_list[id];
        result.push_back(s);
    }
    return result;
}

std::vector<Segment*> get_closest_segment_list_by_LAT_error (std::vector<Segment*> s_list, const double pos[], const double ref_lat, const uint32_t num_size)
{
    std::vector<Segment*> result;
    uint32_t ns = s_list.size();
    double middle_pos[3];
    double min_dist = __DBL_MAX__;
    std::vector<std::pair<double,uint32_t>> error_table;
    for (uint32_t i = 0; i < ns; i++)
    {
        Segment *s = s_list[i];
        s->calc_middle_point(middle_pos);
        double dist = euclidean_norm(middle_pos[0],middle_pos[1],middle_pos[2],pos[0],pos[1],pos[2]);
        double cv = s->calc_propagation_velocity()*M_S_TO_UM_MS;                    // {um/ms}
        double aprox_lat = s->calc_terminal_local_activation_time() + (dist / cv);      // {ms}
        double error = fabs(ref_lat-aprox_lat);
        error_table.push_back(std::make_pair(error,i));
    }
    std::sort(error_table.begin(),error_table.end());
    //for (uint32_t i = 0; i < error_table.size(); i++) // All segments
    // Pick up only the first "N_a" segments
    for (uint32_t i = 0; i < num_size; i++)
    {
        uint32_t id = error_table[i].second;
        Segment *s = s_list[id];
        result.push_back(s);
    }
    return result;
}

Segment* get_closest_segment_by_LAT (std::vector<Segment*> s_list, const double ref_lat)
{
    uint32_t ns = s_list.size();
    double middle_pos[3];
    double min_error = __DBL_MAX__;
    Segment *closest_segment = nullptr;
    for (uint32_t i = 0; i < ns; i++)
    {
        Segment *s = s_list[i];
        double lat = s->lat;
        double error = fabs(ref_lat-lat);

        if ( error < min_error )
        {
            min_error = error;
            closest_segment = s;
        }
    }
    return closest_segment;
}

void print_segment_list (std::vector<Segment*> s_list)
{
    for (uint32_t i = 0; i < s_list.size(); i++)
        s_list[i]->print();
}

/*
void copy_segment_list (std::vector<Node*> in_points, std::vector<Segment*> in, std::vector<Segment*> &out)
{
    for (uint32_t i = 0; i < in.size(); i++)
    {
        uint32_t src_index = in[i]->src->id;
        uint32_t dest_index = in[i]->dest->id;
        double diameter = in[i]->diameter;
        Node *src = in_points[src_index];
        Node *dest = in_points[dest_index];

        Segment *s = new Segment(i,diameter,src,dest,nullptr,nullptr,nullptr);
        s->diameter = in[i]->diameter;
        s->is_mintree = in[i]->is_mintree;

        out.push_back(s);
    }
    // Set the Segment pointers
    for (uint32_t i = 0; i < in.size(); i++)
    {
        Segment *left = in[i]->left;
        Segment *right = in[i]->right;
        Segment *parent = in[i]->parent;

        out[i]->left = (left) ? out[left->id] : nullptr;
        out[i]->right = (right) ? out[right->id] : nullptr;
        out[i]->parent = (parent) ? out[parent->id] : nullptr;
    }
}

void concatenate_segment_list (std::vector<Segment*> in, std::vector<Node*> out_points, std::vector<Segment*> &out,\
                                const uint32_t offset_points, const uint32_t network_id)
{
    const uint32_t offset_segments = out.size();
    for (uint32_t i = 0; i < in.size(); i++)
    {
        uint32_t src_index = in[i]->src->id + offset_points;
        uint32_t dest_index = in[i]->dest->id + offset_points;
        double diameter = in[i]->diameter;
        Node *src = out_points[src_index];
        Node *dest = out_points[dest_index];

        Segment *s = new Segment(i,diameter,src,dest,nullptr,nullptr,nullptr);
        s->id += offset_segments;
        s->diameter = in[i]->diameter;
        s->is_mintree = in[i]->is_mintree;

        out.push_back(s);
    }
    // Update the segment pointers
    for (uint32_t i = 0; i < in.size(); i++)
    {
        Segment *left = in[i]->left;
        Segment *right = in[i]->right;
        Segment *parent = in[i]->parent;

        out[i + offset_segments]->left = (left) ? out[left->id + offset_segments] : nullptr;
        out[i + offset_segments]->right = (right) ? out[right->id + offset_segments] : nullptr;
        out[i + offset_segments]->parent = (parent) ? out[parent->id + offset_segments] : nullptr;
    }
}
*/

Node* generate_bifurcation_node (const uint32_t num_nodes, Segment *iconn)
{
    Node *result = nullptr;
    double middle_pos[3];
    iconn->calc_middle_point(middle_pos);
    
    result = new Node(num_nodes,middle_pos);
    return result;
}

Node* generate_terminal_node (const uint32_t num_nodes, const double pos[])
{
    return (new Node(num_nodes,pos));
}

void move_bifurcation_position (Segment *iconn, Segment *ibiff, Segment *inew, const double pos[])
{
    memcpy(ibiff->dest->pos,pos,sizeof(double)*3);
    memcpy(iconn->src->pos,pos,sizeof(double)*3);
    memcpy(inew->src->pos,pos,sizeof(double)*3);
    iconn->length = iconn->calc_length();
    ibiff->length = ibiff->calc_length();
    inew->length = inew->calc_length();
}

void update_segment_pointers (Segment *iconn, Segment *ibiff, Segment *inew, const bool is_restore)
{
    // Update pointers when building a new segment
    if (!is_restore)
    {
        // iconn:
        iconn->parent = ibiff;
        iconn->src = inew->src;

        // ibiff:
        ibiff->left = inew;     // CONVENTION: Left will always point to terminal
        ibiff->right = iconn;   // CONVENTION: Right will always point to subtree
        if (ibiff->parent != NULL)
        {
            Segment *ibiff_parent = ibiff->parent;
            if (ibiff_parent->left && ibiff_parent->left->id == iconn->id)
                ibiff_parent->left = ibiff;
            if (ibiff_parent->right && ibiff_parent->right->id == iconn->id)
                ibiff_parent->right = ibiff;
        }
    }
    // Update pointers when restoring the tree state
    else
    {
        Segment *ibiff_par = ibiff->parent;

        // iconn:
        iconn->parent = ibiff->parent;
        iconn->src = ibiff->src;

        // ibiff_par:
        if (ibiff_par != NULL)
        {
            if (ibiff_par->right == ibiff)
                ibiff_par->right = iconn;
            if (ibiff_par->left == ibiff)
                ibiff_par->left = iconn;
        }
    }
}

void recalculate_length (std::vector<Segment*> &s_list)
{
    uint32_t ns = s_list.size();
    for (uint32_t i = 0; i < ns; i++)
    {
        Segment *s = s_list[i];
        s->length = s->calc_length();
    }
}

void recalculate_local_activation_time (std::vector<Segment*> &s_list)
{
    uint32_t ns = s_list.size();
    for (uint32_t i = 0; i < ns; i++)
    {
        Segment *s = s_list[i];
        s->lat = s->calc_terminal_local_activation_time();
    }
}

void recalculate_length_branch (Segment *inew, Segment *ibiff)
{
    Segment *tmp = inew;
    while (tmp != ibiff)
    {
        tmp->length = tmp->calc_length();
        tmp = tmp->parent;
    }
}

void recalculate_local_activation_time_branch (Segment *inew, Segment *ibiff)
{
    // Corner case
    if (ibiff->parent)
        ibiff->lat = ibiff->parent->lat + ibiff->calc_local_activation_time();
    
    Segment *tmp = ibiff->left;
    double lat;
    double sum = ibiff->lat;
    while (tmp != inew)
    {
        lat = tmp->calc_local_activation_time();
        tmp->lat = lat + sum;
        sum += lat;
        tmp = tmp->right;
    }
    // Last segment
    lat = tmp->calc_local_activation_time();
    tmp->lat = lat + sum;
}

void compute_distance(std::vector<Segment*> s_list, const double pos[], std::vector< std::pair<double,uint32_t> > &dist_array)
{
    uint32_t ns = s_list.size();
    for (uint32_t i = 0; i < ns; i++)
    {
        double middle_pos[3];
        Segment *s = s_list[i];
        s->calc_middle_point(middle_pos);
        double dist = euclidean_norm(pos[0],pos[1],pos[2],middle_pos[0],middle_pos[1],middle_pos[2]);
        dist_array.push_back(std::make_pair(dist,s->id));
    }
    std::sort(dist_array.begin(),dist_array.end());
}

void write_segment_list (const char filename[], std::vector<Segment*> s_list)
{
    FILE *file = fopen(filename,"w+");
    
    for (uint32_t i = 0; i < s_list.size(); i++)
    {
        Segment *s = s_list[i];

        fprintf(file,"Segment %u (%u,%u) (%g %g %g) - (%g %g %g) -- DIAMETER = %g μm -- LAT = %g ms\n",s->id,s->src->id,s->dest->id,\
                                                                        s->src->pos[0],s->src->pos[1],s->src->pos[2],s->dest->pos[0],s->dest->pos[1],s->dest->pos[2],\
                                                                        s->diameter,s->lat);
        if (s->parent == nullptr)
            fprintf(file,"\tPARENT = NIL\n");
        else
            fprintf(file,"\tPARENT = %u\n",s->parent->id);
        if (s->left == nullptr)
            fprintf(file,"\tLEFT = NIL\n");
        else
            fprintf(file,"\tLEFT = %u\n",s->left->id);
        if (s->right == nullptr)
            fprintf(file,"\tRIGHT = NIL\n");
        else
            fprintf(file,"\tRIGHT = %u\n",s->right->id);
    }

    fclose(file);
}

uint32_t get_total_time_local_activation_terminal ()
{
    return total_time_local_activation_terminal;
}