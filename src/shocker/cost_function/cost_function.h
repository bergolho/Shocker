#ifndef COST_FUNCTION_H
#define COST_FUNCTION_H

#include "../../options/cost_function_config.h"
#include "../network/segment.h"
#include <cassert>

class Evaluation
{
public:
    double eval;
    double term_pos[3];
    bool pass;
    Segment *iconn;
public:
    Evaluation ();
    Evaluation (Segment *iconn, const double eval, const double term_pos[], const bool pass);
    ~Evaluation ();
    void update (Evaluation e);
    void print ();
};

class CostFunction
{
public:
    double beta;
    double alpha;
    double min_degrees_limit;
    double max_degrees_limit;
    double min_segment_length;
    double max_segment_length;
public:
    CostFunction (CostFunctionConfig *cost_function_config);
    void init_parameters (CostFunctionConfig *cost_function_config);
    double eval (std::vector<Segment*> s_list);
    double eval_branch (Segment *inew, Segment *ibiff);
    bool check_evaluations (std::vector<Evaluation> evals, Evaluation &best);
    bool check_restrictions (std::vector<Segment*> s_list, Segment *ibiff);
    bool check_restrictions_pmj (std::vector<Segment*> s_list, Segment *ibiff);
    bool check_angle_restriction (Segment *iconn, Segment *inew);
    bool check_minimum_segment_size (Segment *ibiff);
    bool check_maximum_segment_size (Segment *ibiff);
    bool check_collision (std::vector<Segment*> s_list, Segment *inew);
};

class Segment_Evaluation
{
public:
    double dist;
    double error;
    Segment *iconn;
public:
    Segment_Evaluation (Segment *iconn, const double dist, const double error);
    ~Segment_Evaluation ();
    void print ();
    friend bool operator <(Segment_Evaluation const &a, Segment_Evaluation const &b)
    {
        return a.dist < b.dist;
    }
};

#endif