#include "cost_function.h"

Evaluation::Evaluation ()
{
    this->iconn = nullptr;
    this->eval = __DBL_MAX__;
    this->pass = false;
}

Evaluation::Evaluation (Segment *iconn, const double eval, const double term_pos[], const bool pass)
{
    this->iconn = iconn;
    this->eval = eval;
    memcpy(this->term_pos,term_pos,sizeof(double)*3);
    this->pass = pass;
}

Evaluation::~Evaluation ()
{
    this->iconn = nullptr;
}

void Evaluation::update (Evaluation e)
{
    this->iconn = e.iconn;
    this->eval = e.eval;
    memcpy(this->term_pos,e.term_pos,sizeof(double)*3);
    this->pass = e.pass;
}

void Evaluation::print ()
{
    printf("|| Eval = %g || term_pos = (%g %g %g) || iconn = %u || pass = %d ||\n",this->eval,\
                                                                this->term_pos[0],this->term_pos[1],this->term_pos[2],\
                                                                this->iconn->id,(int)this->pass);
}

CostFunction::CostFunction (CostFunctionConfig *cost_function_config)
{
    init_parameters(cost_function_config);
}

void CostFunction::init_parameters (CostFunctionConfig *cost_function_config)
{
    std::string::size_type sz;
    std::string str;
    cost_function_config->get_parameter_value_from_map("beta",&str);
    beta = std::stod(str,&sz);
    cost_function_config->get_parameter_value_from_map("alpha",&str);
    alpha = std::stod(str,&sz);
    cost_function_config->get_parameter_value_from_map("min_degrees_limit",&str);
    min_degrees_limit = std::stod(str,&sz);
    cost_function_config->get_parameter_value_from_map("max_degrees_limit",&str);
    max_degrees_limit = std::stod(str,&sz);
    cost_function_config->get_parameter_value_from_map("min_segment_length",&str);
    min_segment_length = std::stod(str,&sz);
    cost_function_config->get_parameter_value_from_map("max_segment_length",&str);
    max_segment_length = std::stod(str,&sz);
}

double CostFunction::eval (std::vector<Segment*> s_list)
{
    double eval = 0.0;
    uint32_t ns = s_list.size();
    for (uint32_t i = 0; i < ns; i++)
        eval += s_list[i]->length;
    return eval;
}

double CostFunction::eval_branch (Segment *inew, Segment *ibiff)
{
    Segment *tmp = inew;
    double sum = 0.0;
    while (tmp != ibiff)
    {
        sum += tmp->length;
        tmp = tmp->parent;
    }
    return sum;
}

bool CostFunction::check_restrictions (std::vector<Segment*> s_list, Segment *ibiff)
{
    assert(ibiff->is_bifurcation());

    // ADD RESTRICTION: Check bifurcation size and angle
    Segment *iconn = ibiff->right;
    Segment *inew = ibiff->left;
    //bool has_angle_requirement = check_angle_restriction(iconn,inew);
    
    // ADD RESTRICTION: Check segment sizes
    //bool has_minimum_segment_size = check_minimum_segment_size(ibiff);
    //bool has_maximum_segment_size = check_maximum_segment_size(ibiff);

    // Collision detection:
    bool has_segment_segment_collision = false; 
    while (inew != nullptr)
    {
        has_segment_segment_collision |= check_collision(s_list,inew);
        inew = inew->right;
    }

    //printf("\tangle=%d, min_length=%d, max_length=%d, collision=%d\n",static_cast<int>(has_angle_requirement),\
                                                        static_cast<int>(has_minimum_segment_size),\
                                                        static_cast<int>(has_maximum_segment_size),\
                                                        static_cast<int>(has_segment_segment_collision));

    return !has_segment_segment_collision;
}

bool CostFunction::check_restrictions_pmj (std::vector<Segment*> s_list, Segment *ibiff)
{
    assert(ibiff->is_bifurcation());

    // Check bifurcation size and angle
    Segment *iconn = ibiff->right;
    Segment *inew = ibiff->left;
    //bool has_angle_requirement = check_angle_restriction(iconn,inew);
    
    // Check segment sizes
    //bool has_minimum_segment_size = check_minimum_segment_size(ibiff);
    //bool has_maximum_segment_size = check_maximum_segment_size(ibiff);

    // Collision detection:
    bool has_segment_segment_collision = false; 
    while (inew != nullptr)
    {
        has_segment_segment_collision |= check_collision(s_list,inew);
        inew = inew->right;
    }

    //printf("\tangle=%d, min_length=%d, max_length=%d, collision=%d\n",static_cast<int>(has_angle_requirement),\
                                                        static_cast<int>(has_minimum_segment_size),\
                                                        static_cast<int>(has_maximum_segment_size),\
                                                        static_cast<int>(has_segment_segment_collision));

    return !has_segment_segment_collision;
}

bool CostFunction::check_angle_restriction (Segment *iconn, Segment *inew)
{
    double u[3], v[3];
    
    iconn->calc_unitary_vector(u);
    inew->calc_unitary_vector(v);
    double angle = calc_angle_between_vectors(u,v);

    return (angle > min_degrees_limit && angle < max_degrees_limit);
}

bool CostFunction::check_minimum_segment_size (Segment *ibiff)
{
    double inew_length = 0.0;
    double ibiff_length = 0.0;
    double iconn_length = 0.0;

    // 'inew' branch
    Segment *tmp = ibiff->left;
    while (tmp != nullptr)
    {
        inew_length += tmp->length;
        tmp = tmp->right;
    }

    // 'iconn' branch
    tmp = ibiff->right;
    while (tmp != nullptr && tmp->is_segment())
    {
        iconn_length += tmp->length;
        tmp = tmp->right;
    }

    // 'ibiff' branch
    ibiff_length = ibiff->length;
    tmp = ibiff->parent;
    while (tmp != nullptr && tmp->is_segment())
    {
        ibiff_length += tmp->length;
        tmp = tmp->parent;
    }
    
    return (inew_length > min_segment_length && iconn_length > min_segment_length && ibiff_length > min_segment_length);
}

bool CostFunction::check_maximum_segment_size (Segment *ibiff)
{
    double b_length = 0.0;
    Segment *tmp = ibiff->left;
    while (tmp != nullptr)
    {
        b_length += tmp->length;
        tmp = tmp->right;
    }
    return (b_length < max_segment_length);
}

bool CostFunction::check_collision (std::vector<Segment*> s_list, Segment *inew)
{
    uint32_t ns = s_list.size();

    // Get the reference to the points from the 'inew' segment
    Node *src_inew = inew->src;
    Node *dest_inew = inew->dest;

    bool intersect_1 = false;
    for (uint32_t i = 0; i < ns; i++)
    {
        Segment *s = s_list[i];

        // inew
        if (s->id != inew->id)
            intersect_1 = collision_detection(src_inew->pos[0],src_inew->pos[1],src_inew->pos[2],\
                                            dest_inew->pos[0],dest_inew->pos[1],dest_inew->pos[2],\
                                            s->src->pos[0],s->src->pos[1],s->src->pos[2],\
                                            s->dest->pos[0],s->dest->pos[1],s->dest->pos[2]);
        
        if (intersect_1)
            return true;
    }
    return false;
}

bool CostFunction::check_evaluations (std::vector<Evaluation> evals, Evaluation &best)
{
    uint32_t nevals = evals.size();
    for (uint32_t i = 0; i < nevals; i++)
    {
        bool pass = evals[i].pass; 
        if (pass)
        {
            double eval = evals[i].eval;
            if (eval < best.eval)
            {
                //printf("Best eval = %u\n",i);
                best.update(evals[i]);
                //best.print();
            }
        }
    }
    return (best.iconn) ? true : false;
}

Segment_Evaluation::Segment_Evaluation (Segment *iconn, const double dist, const double error)
{
    this->dist = dist;
    this->error = error;
    this->iconn = iconn;
}

Segment_Evaluation::~Segment_Evaluation ()
{
    this->iconn = nullptr;
}

void Segment_Evaluation::print ()
{
    printf("Dist: %g || Error: %g || Segment: %u\n",this->dist,this->error,this->iconn->id);
}