#include "shocker.h"

struct stop_watch solver_time;
struct stop_watch inactive_points_time;
struct stop_watch active_points_time;
struct stop_watch inactive_search_time;
struct stop_watch active_search_time;
struct stop_watch inactive_evaluation_time;
struct stop_watch active_evaluation_time;
struct stop_watch geodesic_pathway_time;
struct stop_watch check_duplicates_time;
struct stop_watch main_loop_time;
struct stop_watch force_connection_no_lat_error_tolerance_time;
struct stop_watch force_connection_no_distance_criterion_time;
struct stop_watch force_connection_no_distance_criterion_geodesic_time;
struct stop_watch force_connection_no_distance_criterion_line_time;
uint32_t total_time_inactive_points = 0;
uint32_t total_time_active_points = 0;
uint32_t total_time_inactive_search = 0;
uint32_t total_time_active_search = 0;
uint32_t total_time_inactive_evaluation = 0;
uint32_t total_time_active_evaluation = 0;
uint32_t total_time_geodesic_pathway = 0;
uint32_t total_time_check_duplicates = 0;
uint32_t total_time_main_loop = 0;
uint32_t total_time_force_connection_no_lat_error_tolerance = 0;
uint32_t total_time_force_connection_no_distance_criterion = 0;
uint32_t total_time_force_connection_no_distance_criterion_geodesic = 0;
uint32_t total_time_force_connection_no_distance_criterion_line = 0;

Shocker::Shocker(User_Options *options)
{
    this->params = new Config(options);
    this->cloud = new Cloud(options->cloud_config,options->pmj_config,options->output_dir,options->root_pos);
    this->error = new Error();
    this->logger = new Logger(options);
}

Shocker::~Shocker()
{
    if (this->params)
        delete this->params;
    if (this->cloud)
        delete this->cloud;
    if (this->error)
        delete this->error;
    if (this->logger)
        delete this->logger;
    for (uint32_t i = 0; i < this->segment_list.size(); i++)
        delete this->segment_list[i];
    for (uint32_t i = 0; i < this->node_list.size(); i++)
        delete this->node_list[i];
}

void Shocker::grow ()
{
    std::cout << get_color_code(GRAY) << "\n[INFO] Growing Shocker network !" << std::endl;

	init_stop_watch(&solver_time);
    init_stop_watch(&geodesic_pathway_time);    
    
    start_stop_watch(&solver_time);

    // Root placement
    root_placement();

    // Grow the Purkinje network
    grow_tree_using_cloud_points();

    // Pos-processing
    pos_process();    // Prune and write minimum tree     

    long res_time = stop_stop_watch(&solver_time);
	double conv_rate = 1000.0*1000.0*60.0;
    std::cout << get_color_code(GRAY) << "[INFO] Resolution time = " << res_time << " μs (" << res_time/conv_rate << " min)" << std::endl;
    std::cout << get_color_code(GRAY) << "[INFO] Total inactive points time = " << total_time_inactive_points << " μs (" << total_time_inactive_points/conv_rate << " min)" << std::endl;
    std::cout << get_color_code(GRAY) << "[INFO] Total active points time = " << total_time_active_points << " μs (" << total_time_active_points/conv_rate << " min)" << std::endl;
    std::cout << get_color_code(GRAY) << "[INFO] Total inactive points search time = " << total_time_inactive_search << " μs (" << total_time_inactive_search/conv_rate << " min)" << std::endl;
    std::cout << get_color_code(GRAY) << "[INFO] Total active points search time = " << total_time_active_search << " μs (" << total_time_active_search/conv_rate << " min)" << std::endl;
    std::cout << get_color_code(GRAY) << "[INFO] Total inactive evaluation time = " << total_time_inactive_evaluation << " μs (" << total_time_inactive_evaluation/conv_rate << " min)" << std::endl;
    std::cout << get_color_code(GRAY) << "[INFO] Total active evaluation time = " << total_time_active_evaluation << " μs (" << total_time_active_evaluation/conv_rate << " min)" << std::endl;
    std::cout << get_color_code(GRAY) << "[INFO] Total main loop time = " << total_time_main_loop << " μs (" << total_time_main_loop/conv_rate << " min)" << std::endl;
    std::cout << get_color_code(GRAY) << "[INFO] Total force connection without LAT error tolerance time = " << total_time_force_connection_no_lat_error_tolerance << " μs (" << total_time_force_connection_no_lat_error_tolerance/conv_rate << " min)" << std::endl;
    std::cout << get_color_code(GRAY) << "[INFO] Total force connection without distance criterion time = " << total_time_force_connection_no_distance_criterion << " μs (" << total_time_force_connection_no_distance_criterion/conv_rate << " min)" << std::endl;
    std::cout << get_color_code(GRAY) << "[INFO] Total force connection without distance criterion time (geodesic) = " << total_time_force_connection_no_distance_criterion_geodesic << " μs (" << total_time_force_connection_no_distance_criterion_geodesic/conv_rate << " min)" << std::endl;
    std::cout << get_color_code(GRAY) << "[INFO] Total force connection without distance criterion time (line) = " << total_time_force_connection_no_distance_criterion_line << " μs (" << total_time_force_connection_no_distance_criterion_line/conv_rate << " min)" << std::endl;
    std::cout << get_color_code(GRAY) << "[INFO] Total geodesic pathway time = " << total_time_geodesic_pathway << " μs (" << total_time_geodesic_pathway/conv_rate << " min)" << std::endl;
    std::cout << get_color_code(GRAY) << "[INFO] Total check duplicates time = " << total_time_check_duplicates << " μs (" << total_time_check_duplicates/conv_rate << " min)" << std::endl;
    std::cout << get_color_code(GRAY) << "[INFO] Total local activation terminal time = " << get_total_time_local_activation_terminal() << " μs (" << get_total_time_local_activation_terminal()/conv_rate << " min)" << std::endl;

    // Output network information
    print_network_info();
    
    // Write information about the simulation to file
    logger->write_simulation_time(this->params->output_dir,res_time,total_time_inactive_points,total_time_active_points,\
                                total_time_inactive_search,total_time_active_search,\
                                total_time_inactive_evaluation,total_time_active_evaluation,\
                                total_time_main_loop,total_time_force_connection_no_lat_error_tolerance,total_time_force_connection_no_distance_criterion,\
                                total_time_force_connection_no_distance_criterion_geodesic,total_time_force_connection_no_distance_criterion_line,\
                                total_time_geodesic_pathway,\
                                conv_rate);
}

void Shocker::grow_tree_using_cloud_points ()
{
    bool sucess = false;
    bool using_pmj_location = this->params->using_pmj_location;
    uint32_t connection_rate = this->cloud->pmj_data->connection_rate;

    // Get reference to the pointers
    CostFunction *cost_fn = this->params->cost_fn;
    Logger *logger = this->logger;

    start_stop_watch(&main_loop_time);

    // ===============================================================================================================
    // MAIN ITERATION LOOP
    //while (this->params->num_terminals < this->params->N_term)
    while (this->cloud->cloud_data->counter_cloud_passes == 0)
    {
        std::cout << get_color_code(GRAY) << PRINT_LINE << std::endl;
        std::cout << get_color_code(GRAY) << "[shocker] Working on terminal number " << this->params->num_terminals+1 << std::endl;

        // Generate a new terminal point in the network
        if (!using_pmj_location)
            sucess = generate_terminal();
        else
        {
            sucess = generate_terminal();

            // Process the PMJs if the current number of terminals is multiple of the connection_rate
            if (this->params->num_terminals % connection_rate == 0)
            {
                sucess = attempt_pmj_connection();
            }   
        }  

        std::cout << get_color_code(GRAY) << PRINT_LINE << std::endl;
    }
    // ===============================================================================================================
    // Write the full network when the main loop is finished
    //write_full_network_to_vtk(); 

    total_time_main_loop += stop_stop_watch(&main_loop_time);                      
}

void Shocker::root_placement ()
{
    bool using_initial_network = this->params->using_initial_network;
    (using_initial_network) ? load_network_state() : make_root_default();
    update_and_check_integrity();
    logger->write_full_network_to_vtk(this->node_list,this->segment_list,this->params->output_dir,this->params->his_offset,this->params->num_terminals);
}

void Shocker::make_root_default ()
{
    std::cout << get_color_code(GRAY) << "[INFO] Making root using default method" << std::endl;
    double x_prox[3], x_dist[3];
    Cloud_Data *cloud_data = this->cloud->cloud_data;
    double l_d = this->params->l_d;
    memcpy(x_prox,this->params->root_pos,sizeof(double)*3);

    cloud_data->calc_distal_root_position(x_prox,l_d,x_dist);
    generate_root_segment(x_prox,x_dist);
}

void Shocker::generate_root_segment (const double x_prox[], const double x_dist[])
{
    Surface_Data *surface_data = this->cloud->surface_data;
    std::vector<CloudPoint> geodesic_points = surface_data->compute_geodesic_pathway(x_prox,x_dist,true,nullptr);
    
    for (uint32_t i = 0; i < geodesic_points.size(); i++)
    {
        double *node_pos = geodesic_points[i].pos;
        Node *node = new Node(i,node_pos);
        this->node_list.push_back(node);
    }
    double d = this->params->start_diameter;
    for (uint32_t i = 0; i < geodesic_points.size()-1; i++)
    {
        Segment *segment = new Segment(i,d,0,false,this->node_list[i],this->node_list[i+1],nullptr,nullptr,nullptr);
        this->segment_list.push_back(segment);
    }
    for (uint32_t i = 0; i < this->segment_list.size()-1; i++)
    {
        this->segment_list[i]->right = this->segment_list[i+1];
        this->segment_list[i+1]->parent = this->segment_list[i];
    }

    this->params->num_terminals = 1;
    recalculate_length(this->segment_list);
    recalculate_local_activation_time(this->segment_list);
}

bool Shocker::generate_terminal ()
{
    uint32_t tosses = 0;
    bool sucess = false;
    bool point_is_ok = false;
    double his_offset = this->params->his_offset;
    double l_min = calc_lmin(this->params->l_d,this->params->num_terminals);
    Cloud_Data *cloud_data = this->cloud->cloud_data;
    CloudPoint *p = nullptr;
    Evaluation best;
    double branch_size;
    
    start_stop_watch(&inactive_points_time);

    while (!point_is_ok)
    {
        // Array of feasible segments to connect the new terminal
        std::vector<Segment*> feasible_segments;

        // Sort a terminal position from the region cloud
        p = cloud_data->sort_point();
        std::cout << get_color_code(GRAY) << "[INFO] Selected index " << p->id << " = (" << p->pos[0] << " " << p->pos[1] << " " << p->pos[2] << ")" << std::endl;

        start_stop_watch(&inactive_search_time);
        
        // Fill the feasible segments array
        bool cs_flag = connection_search(p->pos,l_min);
        bool ffs_flag = fill_feasible_segments(feasible_segments,p->pos);  
        
        total_time_inactive_search += stop_stop_watch(&inactive_search_time);
        
        // The cloud point must fulfil the distance criterion and have at least one feasible segment to connect
        if (cs_flag && ffs_flag)
        {
            start_stop_watch(&inactive_evaluation_time);

            // Evaluate all possible configurations with the cost function
            point_is_ok = evaluate_cost_function_cloud_point(feasible_segments,p,best);
            
            total_time_inactive_evaluation += stop_stop_watch(&inactive_evaluation_time);
        }

        // If there is no feasible segment to connect the terminal point
        // we update minimum distance and sort another point
        if (!point_is_ok)
            l_min = update_minimum_distance(tosses,l_min);

        // Error handling
        if (check_minimum_distance(l_min))
        {
            fprintf(stderr,"[INFO] ERROR! Minimum distance value below acceptable value reached!\n");
            exit(EXIT_FAILURE);
        } 
    }

    // Build the new segment at the best position found in the evaluations
    Segment *inew = build_branch(best.iconn,p->pos,branch_size);
    
    // Check for errors
    update_and_check_integrity();

    total_time_inactive_points += stop_stop_watch(&inactive_points_time);

    // Write the current stage of the Purkinje network to a file
    logger->write_full_network_to_vtk(this->node_list,this->segment_list,this->params->output_dir,this->params->his_offset,this->params->num_terminals);

    return point_is_ok;
}

bool Shocker::generate_pmj (PMJPoint *pmj)
{
    uint32_t tosses = 0;
    bool sucess = false;
    bool pmj_is_ok = false;
    double his_offset = this->params->his_offset;
    double l_min = calc_lmin(this->params->l_d,this->params->num_terminals);
    Cloud_Data *cloud_data = this->cloud->cloud_data;
    Evaluation best;
    double branch_size;

    // Array of feasible segments to connect the PMJ
    std::vector<Segment*> feasible_segments;

    start_stop_watch(&active_points_time);
    start_stop_watch(&active_search_time);

    // Fill the feasible segments array
    bool cs_flag = connection_search(pmj->pos,l_min);
    //bool ffs_flag = fill_feasible_segments(feasible_segments,pmj->pos);
    bool ffs_flag = fill_feasible_segments_pmj(feasible_segments,pmj->pos,pmj->ref_value);
    
    total_time_active_search += stop_stop_watch(&active_search_time);
    
    // The PMJ point must have at least one feasible segment to connect
    if (cs_flag && ffs_flag)
    {
        start_stop_watch(&active_evaluation_time);

        // Evaluate all possible configurations with the cost function
        pmj_is_ok = evaluate_cost_function_pmj(feasible_segments,pmj,best);

        total_time_active_evaluation += stop_stop_watch(&active_evaluation_time);
    }
    else
    {
        std::cerr << get_color_code(RED) << "\t[INFO] ERROR! No feasible segment found for PMJ point " << pmj->id << std::endl;
        if (!cs_flag)
            std::cerr << get_color_code(RED) << "\t[INFO] ERROR! Distance criterion not attended for PMJ point " << pmj->id << std::endl;
        if (!ffs_flag)
            std::cerr << get_color_code(RED) << "\t[INFO] ERROR! Fill feasible segments returned no segments for PMJ point " << pmj->id << std::endl;
    }

    // If there is no feasible segment or valid configuration to connect the PMJ point
    // we return a FAILURE signal
    if (!pmj_is_ok)
    {
        if (!best.iconn && cs_flag && ffs_flag)
            std::cerr << get_color_code(RED) << "\t[INFO] ERROR! All feasible segments generate collision for PMJ point " << pmj->id << std::endl;
        std::cerr << get_color_code(RED) << "\t[INFO] ERROR! Restriction failure for PMJ point " << pmj->id << std::endl;
        total_time_active_points += stop_stop_watch(&active_points_time);
        return pmj_is_ok;
    } 
    // Otherwise, we connect the PMJ to tree and check the LAT error tolerance
    else
    {
        // Build the new segment at the best segment found in the evaluations
        Segment *inew = build_branch(best.iconn,pmj->pos,branch_size);
        
        // Check if the PMJ is connected within the LAT error tolerance
        pmj_is_ok = evaluate_pmj_local_activation_time(inew,pmj);
        if (pmj_is_ok) 
        {
            std::cout << get_color_code(GREEN) << "\t[INFO] PMJ point " << pmj->id << " was connected!" << std::endl;
            inew->tag_pathway();
            update_and_check_integrity();
            if (this->cloud->pmj_data->lat_error_tolerance != __DBL_MAX__)
                this->error->counter_main_loop_connections++;
            else
                this->error->counter_no_lat_error_tolerance_connections++;

            // Write the current stage of the Purkinje network to a file
            logger->write_full_network_to_vtk(this->node_list,this->segment_list,this->params->output_dir,this->params->his_offset,this->params->num_terminals);
        }
        else
        {
            std::cerr << get_color_code(RED) << "\t[INFO] ERROR! LAT error tolerance was not met for PMJ point " << pmj->id << std::endl;
        }
        total_time_active_points += stop_stop_watch(&active_points_time);
        return pmj_is_ok;
    }     
}

bool Shocker::attempt_pmj_connection ()
{
    uint32_t cur_num_connected_pmjs = 0;
    uint32_t prev_num_connected_pmjs = 0;
    uint32_t total_num_pmjs_connected_before = this->cloud->pmj_data->total_num_connected;

    // Process all the PMJs until we receive a FAILURE signal in everyone that is not connected
    bool sucess;
    do
    {
        sucess = process_pmjs(cur_num_connected_pmjs,prev_num_connected_pmjs);
                        
    } while (prev_num_connected_pmjs != cur_num_connected_pmjs);

    return (this->cloud->pmj_data->total_num_connected > total_num_pmjs_connected_before) ? true : false;
}

bool Shocker::process_pmjs (uint32_t &prev_num_connected_pmjs, uint32_t &cur_num_connected_pmjs)
{
    bool sucess = false;
    PMJ_Data *pmj_data = this->cloud->pmj_data;
    prev_num_connected_pmjs = cur_num_connected_pmjs;

    std::cout << get_color_code(GRAY) << PRINT_DOTS << std::endl;
    for (uint32_t i = 0; i < pmj_data->pmjs.size(); i++)
    {
        PMJPoint *pmj = pmj_data->pmjs[i];
        if (!pmj->connected)
        {
            std::cout << get_color_code(GRAY) << "[shocker] Trying to connect PMJ point " << pmj->id << " ..." << std::endl;
            bool ret = generate_pmj(pmj);
            pmj->connected = ret;

            if (ret)
            {
                std::cout << get_color_code(GREEN) << "\t\t\t[!] SUCESS!" << std::endl;
                pmj_data->total_num_connected++;
                cur_num_connected_pmjs++;
            }
            sucess |= ret;
        }
    }
    std::cout << get_color_code(GRAY) << PRINT_DOTS << std::endl;
    return sucess;
}

bool Shocker::connection_search (const double pos[], const double l_min)
{
    double region_radius = this->cloud->cloud_data->REGION_RADIUS;
    for (uint32_t i = 0; i < this->segment_list.size(); i++)
    {
        Segment *s = this->segment_list[i];
        if (!distance_criterion(s,pos,l_min))
            return false;
    }
    return true;
}

bool Shocker::fill_feasible_segments (std::vector<Segment*> &feasible_segments, const double pos[])
{
    uint32_t ns = this->segment_list.size();
    uint32_t N_p = this->params->np;
    std::vector< std::pair<double,uint32_t> > arr;
    // TODO: Use a Kd-tree to calculate the distance to the closest point using a range subroutine with Ncon as a parameter
    // FunctionName: vtkKdtree::FindClosestNPoints()
    for (uint32_t i = 0; i < ns; i++)
    {
        Segment *s = this->segment_list[i];
        double middle_pos[3];
        s->calc_middle_point(middle_pos);
        double value = euclidean_norm(middle_pos[0],middle_pos[1],middle_pos[2],pos[0],pos[1],pos[2]);
        arr.push_back(std::make_pair(value,s->id));
    }
    // Sort the feasible segments by their distance
    std::sort(arr.begin(),arr.end());
    // Fill the 'feasible_segments' array with N_p segments
    for (uint32_t i = 0; i < arr.size() && feasible_segments.size() < N_p; i++)
    {
        uint32_t id = arr[i].second;
        Segment *s = this->segment_list[id];
        feasible_segments.push_back(s);
    }
    return (feasible_segments.size() > 0) ? true : false;
}

bool Shocker::fill_feasible_segments_pmj (std::vector<Segment*> &feasible_segments, const double pos[], const double ref_lat)
{
    uint32_t ns = this->segment_list.size();
    uint32_t N_a = this->params->na;
    std::vector< std::pair<double,uint32_t> > arr;
    double his_offset = this->params->his_offset;
    // TODO: Use a Kd-tree to calculate the distance to the closest point using a range subroutine with Ncon as a parameter
    // FunctionName: vtkKdtree::FindClosestPoint()
    for (uint32_t i = 0; i < ns; i++)
    {
        Segment *s = this->segment_list[i];
        double middle_pos[3];
        s->calc_middle_point(middle_pos);
        double dist = euclidean_norm(middle_pos[0],middle_pos[1],middle_pos[2],pos[0],pos[1],pos[2]);  
        double cv = s->calc_propagation_velocity()*M_S_TO_UM_MS;                    // {um/ms}
        double aprox_lat = s->calc_terminal_local_activation_time() + (dist / cv);      // {ms}
        double error = fabs(ref_lat-aprox_lat);
        arr.push_back(std::make_pair(error,s->id));
    }
    // Sort the feasible segments by their distance
    std::sort(arr.begin(),arr.end());
    for (uint32_t i = 0; i < arr.size() && feasible_segments.size() < N_a; i++)
    {
        uint32_t id = arr[i].second;
        Segment *s = this->segment_list[id];
        feasible_segments.push_back(s);
    }
    return (feasible_segments.size() > 0) ? true : false;
}

bool Shocker::evaluate_cost_function_cloud_point (std::vector<Segment*> &feasible_segments, CloudPoint *point, Evaluation &best)
{
    CostFunction *cost_fn = this->params->cost_fn;
    std::vector<Evaluation> evaluations;
    double *term_pos = point->pos;
    uint32_t num_fs = feasible_segments.size();
    //uint32_t num_fsi = this->params->pp*this->params->nconn;
    //uint32_t num_size = (num_fs < num_fsi) ? num_fs : num_fsi;
    double branch_size;

    // We only use a percentage of the NCONN segments
    for (uint32_t i = 0; i < num_fs; i++)
    {
        //printf("\tEvaluating segment %u [%u/%u]\n",feasible_segments[i]->id,i,num_fs);

        // Build a new branch with a bifurcation point at the middle position of 'iconn'
        Segment *iconn = feasible_segments[i];
        Segment *inew = build_branch(iconn,point->pos,branch_size);
        Segment *ibiff = iconn->parent;

        if (inew)
        {
            // COST FUNCTION: Minimize network length (Optimized)
            // Evaluate the cost function only for the current branch
            //double eval_branch = cost_fn->eval_branch(inew, ibiff);
            //double eval = base_eval + eval_branch;
            double eval = branch_size;
            bool pass = cost_fn->check_restrictions(this->segment_list,ibiff);

            // COST FUNCTION: Minimize network length
            // Evaluate the cost function
            //double eval = cost_fn->eval(this->segment_list);
            //bool pass = cost_fn->check_restrictions(this->segment_list,ibiff);

            Evaluation e(iconn,eval,term_pos,pass);
            evaluations.push_back(e);

            // Restore the network to its state before the new terminal
            prune_branch(inew);
        }
    }

    // Return TRUE if we have at least one feasible segment to make the connection
    return cost_fn->check_evaluations(evaluations,best) ? true : false;
}

bool Shocker::evaluate_cost_function_pmj (std::vector<Segment*> &feasible_segments, PMJPoint *pmj, Evaluation &best)
{
    CostFunction *cost_fn = this->params->cost_fn;
    std::vector<Evaluation> evaluations;
    double *term_pos = pmj->pos;
    double his_offset = this->params->his_offset;
    double branch_size;
    uint32_t num_fs = feasible_segments.size();
    //uint32_t num_fsa = this->params->pa*this->params->nconn;
    //uint32_t num_size = (num_fs < num_fsa) ? num_fs : num_fsa;

    // We only use a percentage of the NCONN segments
    for (uint32_t i = 0; i < num_fs; i++)
    {
        // Build a new branch with a bifurcation point at the middle position of 'iconn'
        Segment *iconn = feasible_segments[i];
        Segment *inew = build_branch(iconn,pmj->pos,branch_size);
        Segment *ibiff = iconn->parent;

        if (inew)
        {
            // COST FUNCTION: Minimize LAT error
            // Evaluate the cost function
            double ref_lat = pmj->ref_value;
            double aprox_lat = inew->calc_terminal_local_activation_time() + his_offset;
            //double aprox_lat = inew->lat + his_offset;
            double eval = fabs(ref_lat-aprox_lat);
            bool pass = cost_fn->check_restrictions_pmj(this->segment_list,ibiff);

            Evaluation e(iconn,eval,term_pos,pass);
            evaluations.push_back(e);

            // Restore the network to its state before the new terminal
            prune_branch(inew);
        }
    }
    // Return TRUE if we have at least one feasible segment to make the connection
    return cost_fn->check_evaluations(evaluations,best) ? true : false;
}

bool Shocker::evaluate_pmj_local_activation_time (Segment *inew, PMJPoint *pmj)
{
    PMJ_Data *pmj_data = this->cloud->pmj_data;
    double his_offset = this->params->his_offset;
    double lat_error_tolerance = this->cloud->pmj_data->lat_error_tolerance;

    // Calculate the LAT error for the PMJ
    double ref_value = pmj->ref_value;
    double aprox_value = inew->calc_terminal_local_activation_time() + his_offset;
    //double aprox_value = inew->lat + his_offset;
    double lat_error = (ref_value - aprox_value);
    std::cout << get_color_code(GRAY) << "\t[PMJ " << pmj->id << "] Ref: " << ref_value << " ms || Aprox: " << aprox_value << "ms || Error: " << lat_error << " ms" << std::endl;

    // Check if the LAT error is less then the absolute error tolerance
    if (fabs(lat_error) < lat_error_tolerance)
    {
        pmj->aprox_value = aprox_value;
        pmj->error = lat_error;
        inew->dest->is_active = true;
        return true;
    }
    // If the error is above the tolerance we prune the segment 
    else
    {
        prune_branch(inew);
        return false;
    }
}

bool Shocker::distance_criterion (Segment *s, const double pos[], const double l_min)
{
    double d_proj = s->calc_dproj(pos);

    double d_crit;
    if (d_proj >= 0 && d_proj <= 1)
        d_crit = s->calc_dortho(pos);
    else
        d_crit = s->calc_dend(pos);

    return (d_crit < l_min) ? false : true;
}

bool Shocker::prune_inactive_segments ()
{
    // Tag the segments of the minimum network
    uint32_t num_tagged_segments = tag_minimum_network_segments();

    // Save the nodes and segments that are part of the minimum network
    std::map<uint32_t,uint32_t> nodes_map;
    std::map<uint32_t,uint32_t> segments_map;
    build_minimum_network_maps(nodes_map,segments_map);

    std::vector<Node> min_network_nodes;
    std::vector<Segment> min_network_segments;
    build_minimum_network_topology(nodes_map,segments_map,min_network_nodes,min_network_segments);

    // Reset the network to its initial state
    reset_network_state();

    // Copy the points and segments of the minimum network
    for (uint32_t i = 0; i < min_network_nodes.size(); i++)
    {
        double pos[3];
        uint32_t id = min_network_nodes[i].id;
        pos[0] = min_network_nodes[i].pos[0];
        pos[1] = min_network_nodes[i].pos[1];
        pos[2] = min_network_nodes[i].pos[2];
        bool is_active = min_network_nodes[i].is_active;
        
        Node *n = new Node(id,pos,is_active);
        this->node_list.push_back(n);
    }
    for (uint32_t i = 0; i < min_network_segments.size(); i++)
    {
        uint32_t id = min_network_segments[i].id;
        double diameter = min_network_segments[i].diameter;
        double lat = min_network_segments[i].lat;
        bool min_tree = min_network_segments[i].mintree;
        uint32_t src_id = min_network_segments[i].src->id;
        uint32_t dest_id = min_network_segments[i].dest->id;
        Node *src = this->node_list[src_id];
        Node *dest = this->node_list[dest_id];

        Segment *s = new Segment(id,diameter,lat,min_tree,src,dest,NULL,NULL,NULL);
        this->segment_list.push_back(s);
    }
    for (uint32_t i = 0; i < min_network_segments.size(); i++)
    {
        Segment *parent = min_network_segments[i].parent;
        Segment *left = min_network_segments[i].left;
        Segment *right = min_network_segments[i].right;
        int par_id = (parent) ? parent->id : -1;
        int left_id = (left) ? left->id : -1;
        int right_id = (right) ? right->id : -1;

        this->segment_list[i]->parent = (par_id != -1) ? this->segment_list[par_id] : NULL;
        this->segment_list[i]->left = (left_id != -1) ? this->segment_list[left_id] : NULL;
        this->segment_list[i]->right = (right_id != -1) ? this->segment_list[right_id] : NULL;
    }
    this->params->num_terminals = this->cloud->pmj_data->total_num_connected;

    recalculate_length(this->segment_list);
    recalculate_local_activation_time(this->segment_list);

    return true;
}

bool Shocker::check_minimum_distance (const double l_min)
{
    return (l_min < L_MIN_LIMIT) ? true : false;
}

bool Shocker::pos_process ()
{
    bool sucess = false;
    bool using_pmj_location = this->params->using_pmj_location;
    PMJ_Data *pmj_data = this->cloud->pmj_data;
    uint32_t total_num_pmjs = pmj_data->pmjs.size();
    Logger *logger = this->logger;
    
    // Connect any remaining PMJs to the closest segment in the network
    if (using_pmj_location)
    {
        // If not a single PMJ is connected, prune the entire tree and call the grow function again!
        while (pmj_data->total_num_connected == 0)
        {
            std::cout << get_color_code(PURPLE) << "[!] No PMJ was connected in this run!" << std::endl;
            prune_inactive_segments();
            
            bool using_initial_network = this->params->using_initial_network;
            (using_initial_network) ? load_network_state() : make_root_default();
            update_and_check_integrity();
            //logger->write_full_network_to_vtk(this->node_list,this->segment_list,this->params->output_dir,this->params->his_offset,this->params->num_terminals);

            // Reset the cloud counter passes
            this->cloud->cloud_data->counter_cloud_passes = 0;
            this->cloud->cloud_data->cur_id = 0;
            
            grow_tree_using_cloud_points();
        }

        // Reduce to the minimum network
        prune_inactive_segments();

        // Check if all PMJs are connected within the LAT error tolerance
        if (pmj_data->total_num_connected != total_num_pmjs)
        {
            start_stop_watch(&force_connection_no_lat_error_tolerance_time);
            
            // Drop the LAT error tolerance and try to connect the remaining PMJs to the segment with minimum error
            // REMAINDER: the distance criterion is still active in this function!
            sucess = force_connection_without_lat_error_tolerance();

            total_time_force_connection_no_lat_error_tolerance += stop_stop_watch(&force_connection_no_lat_error_tolerance_time);

            if (sucess)
            {
                std::cout << get_color_code(GREEN) << "[!] All PMJs are connected to the network!" << std::endl;
            }
            else
            {
                start_stop_watch(&force_connection_no_distance_criterion_time);

                // Now, drop the distance criterion
                sucess = force_connection_without_distance_criterion();

                total_time_force_connection_no_distance_criterion += stop_stop_watch(&force_connection_no_distance_criterion_time);

                if (sucess)
                {
                    std::cout << get_color_code(GREEN) << "[!] All PMJs are connected to the network!" << std::endl;
                }
                else
                {
                    std::cout << get_color_code(RED) << "[-] ERROR! Not all PMJs are connected to the network!" << std::endl;
                    exit(EXIT_FAILURE);
                }
            }
        }
        write_minimum_network_to_vtk();
    }
    return true;
}

bool Shocker::force_connection_without_distance_criterion ()
{
    std::cout << get_color_code(YELLOW) << "[shocker] Force connection of the remaining PMJs to the best segment in the network without the distance criterion" << std::endl;

    PMJ_Data *pmj_data = this->cloud->pmj_data;
    CostFunction *cost_fn = this->params->cost_fn;
    uint32_t N_a = this->params->na;
    double his_offset = this->params->his_offset;
    double branch_size;
    for (uint32_t i = 0; i < pmj_data->pmjs.size(); i++)
    {
        PMJPoint *pmj = pmj_data->pmjs[i];
        Evaluation best;
        if (!pmj->connected)
        {
            // TODO: Build the two closest "segment_list" together in the same loop
            // Sort all the segments by the distance and LAT error to the unconnected PMJ point
            // Return only the first "N_a" segments in each array
            std::vector<Segment*> nearest_list_by_dist = get_closest_segment_list_by_distance(this->segment_list,pmj->pos, N_a);
            std::vector<Segment*> nearest_list_by_LAT = get_closest_segment_list_by_LAT_error(this->segment_list,pmj->pos,pmj->ref_value, N_a);
            uint32_t num_nl = nearest_list_by_LAT.size();

            // Evaluate the active cost function on the first "N_a" segments of the tree 
            bool pmj_is_ok = evaluate_cost_function_pmj(nearest_list_by_LAT,pmj,best);

            // If we find a suitable segment to build a geodesic pathway we build the branch
            // The best segment will be stored in the "best" variable
            if (pmj_is_ok)
            {
                start_stop_watch(&force_connection_no_distance_criterion_geodesic_time);

                // Build the new segment at the best segment found in the evaluations
                Segment *inew = build_branch(best.iconn,pmj->pos,branch_size);
                
                // Check if the PMJ is connected within the LAT error tolerance
                pmj_is_ok = evaluate_pmj_local_activation_time(inew,pmj);

                // If the point is valid, we tag the pathway until the root as part of the minimum network
                if (pmj_is_ok) 
                {
                    inew->tag_pathway();
                    update_and_check_integrity();
                    this->error->counter_no_distance_criterion_connections++;
                    this->error->counter_no_distance_criterion_connections_geodesic++;

                    //printf("\t[DEBUG] The real LAT error of PMJ %u was equal to %g ms. Segment id which the PMJ was connected is %u\n",pmj->id,lat_error,nearest_list_by_LAT[j]->id);

                    std::cout << get_color_code(YELLOW) << "\t[PMJ " << pmj->id << "] Ref: " << pmj->ref_value << " ms || Aprox: " << pmj->aprox_value << "ms || Error: " << pmj->error << " ms" << std::endl;
                    std::cout << get_color_code(PURPLE) << "\tConnected with geodesic pathway!" << std::endl;
                    std::cout << get_color_code(GREEN) << "\t\t\t[!] SUCESS!" << std::endl;
                    pmj_data->total_num_connected++;
                }   

                total_time_force_connection_no_distance_criterion_geodesic += stop_stop_watch(&force_connection_no_distance_criterion_geodesic_time);
            }
            // Otherwise, we could not generate a geodesic pathway to the PMJ. This could happen because of proximity or collision
            // In this case, we will search for the closest segment by distance to the PMJ a select the one with minimum LAT error
            // This connection will be made using a straight line
            else
            {
                start_stop_watch(&force_connection_no_distance_criterion_line_time);

                // Compute the LAT error for the first "N_a" closest segments to the PMJ
                // Since the final segment is a straight line, the LAT error will be the real one in this case
                const uint32_t nclosest = 5;
                std::vector<std::pair<double,uint32_t>> arr;
                //for (uint32_t j = 0; j < nearest_list_by_dist.size(); j++)
                for (uint32_t j = 0; j < nclosest; j++)
                {
                    Segment *s = nearest_list_by_dist[j];
                    double middle_pos[3];
                    s->calc_middle_point(middle_pos);
                    double dist = euclidean_norm(middle_pos[0],middle_pos[1],middle_pos[2],pmj->pos[0],pmj->pos[1],pmj->pos[2]);  
                    double cv = s->calc_propagation_velocity()*M_S_TO_UM_MS;                        // {um/ms}
                    double aprox_lat = s->calc_terminal_local_activation_time() + (dist / cv);      // {ms}
                    double error = fabs(pmj->ref_value-aprox_lat);
                    arr.push_back(std::make_pair(error,s->id));
                }   

                // We only consider a certain number of closest segments in order to minimize having straight lines out of the endocardium surface
                // Get from the first NCLOSEST segments, the one that returns the minimum LAT error
                Segment *nearest = nullptr;
                //uint32_t nclosest = nearest_list_by_dist.size();
                double min_error = __DBL_MAX__;
                for (uint32_t j = 0; j < nclosest; j++)
                {
                    double error = arr[j].first;
                    if (error < min_error)
                    {
                        min_error = error;
                        nearest = nearest_list_by_dist[j];
                    }
                }

                // Make the necessary adjustments to include the new line to the network
                nearest->calc_middle_point(nearest->middle_pos);

                uint32_t nn = this->node_list.size();
                uint32_t ns = this->segment_list.size();
                double d = nearest->diameter;

                double *node_pos = nearest->middle_pos;
                double *pmj_pos = pmj->pos;
                Node *biff_node = new Node(nn,node_pos);
                Node *pmj_node = new Node(nn+1,pmj_pos);
                this->node_list.push_back(biff_node);
                this->node_list.push_back(pmj_node);

                Segment *ibiff = new Segment(ns,d,0,false,nearest->src,biff_node,nullptr,nullptr,nearest->parent);
                Segment *inew = new Segment(ns+1,d,0,false,biff_node,pmj_node,nullptr,nullptr,nearest);
                this->segment_list.push_back(ibiff);
                this->segment_list.push_back(inew);

                nearest->parent = ibiff;
                nearest->src = biff_node;
                ibiff->right = nearest;
                ibiff->left = inew;

                // Calculate the LAT error for the PMJ
                double ref_value = pmj->ref_value;
                double aprox_value = inew->calc_terminal_local_activation_time() + his_offset;
                //double aprox_value = inew->lat + his_offset;
                double lat_error = (ref_value - aprox_value);

                pmj->aprox_value = aprox_value;
                pmj->error = lat_error;
                pmj->connected = true;
                inew->dest->is_active = true;

                // Update the network parameters
                this->params->num_terminals++;
                recalculate_length(this->segment_list);
                recalculate_local_activation_time(this->segment_list);

                inew->tag_pathway();
                update_and_check_integrity();
                this->error->counter_no_distance_criterion_connections++;
                this->error->counter_no_distance_criterion_connections_line++;

                std::cout << get_color_code(YELLOW) << "\t[PMJ " << pmj->id << "] Ref: " << ref_value << " ms || Aprox: " << aprox_value << "ms || Error: " << lat_error << " ms" << std::endl;
                std::cout << get_color_code(PURPLE) << "\tConnected with straight line!" << std::endl;
                std::cout << get_color_code(GREEN) << "\t\t\t[!] SUCESS!" << std::endl;
                pmj_data->total_num_connected++;

                total_time_force_connection_no_distance_criterion_line += stop_stop_watch(&force_connection_no_distance_criterion_line_time);
            }
            logger->write_full_network_to_vtk(this->node_list,this->segment_list,this->params->output_dir,this->params->his_offset,this->params->num_terminals);
        }
    }
    return (pmj_data->total_num_connected == pmj_data->pmjs.size());
}

bool Shocker::force_connection_without_lat_error_tolerance ()
{
    std::cout << get_color_code(YELLOW) << "[shocker] Force connection of the remaining PMJs to the best segment in the network without the LAT error tolerance" << std::endl;

    PMJ_Data *pmj_data = this->cloud->pmj_data;
    double his_offset = this->params->his_offset;
    for (uint32_t i = 0; i < pmj_data->pmjs.size(); i++)
    {
        PMJPoint *pmj = pmj_data->pmjs[i];
        if (!pmj->connected)
        {
            // Eliminate the LAT error tolerance
            // REMAINDER: Here we still considered the first "N_a" segments
            pmj_data->lat_error_tolerance = __DBL_MAX__;

            // Increase "N_a" to be the total number of segments
            //this->params->na = this->segment_list.size();
            
            std::cout << get_color_code(YELLOW) << "[shocker] Trying to connect PMJ point " << pmj->id << " ..." << std::endl;
            bool ret = generate_pmj(pmj);
            pmj->connected = ret;

            if (ret)
            {
                std::cout << get_color_code(GREEN) << "\t\t\t[!] SUCESS!" << std::endl;
                pmj_data->total_num_connected++;
            }
        }
    }
    return (pmj_data->total_num_connected == pmj_data->pmjs.size());
}

uint32_t Shocker::tag_minimum_network_segments ()
{
    PMJ_Data *pmj_data = this->cloud->pmj_data;
    uint32_t num_pmjs = pmj_data->pmjs.size();
    uint32_t ns = this->segment_list.size();

    for (uint32_t i = 0 ; i < ns; i++)
    {
        Segment *s = this->segment_list[i];
        Node *dest = s->dest;
        if (dest->is_active)
        {
            Segment *tmp = s;
            while (tmp != nullptr)
            {
                // Mark all the segments that are in the pathway to the root
                tmp->mintree = true;
                tmp = tmp->parent;
            }
        }
    }

    // Count the total number of tagged segments
    uint32_t num_tagged_segments = 0;
    for (uint32_t i = 0 ; i < ns; i++)
    {
        Segment *s = this->segment_list[i];
        if (s->is_mintree()) num_tagged_segments++;
    }

    return num_tagged_segments;
}

Segment* Shocker::search_active_terminal_segment (PMJPoint *pmj)
{
    Segment *result = nullptr;
    double min_dist = __DBL_MAX__;
    uint32_t ns = this->segment_list.size();
    for (uint32_t i = 0 ; i < ns; i++)
    {
        Segment *s = this->segment_list[i];
        Node *x_distal = s->dest;
        double dist = euclidean_norm(x_distal->pos[0],x_distal->pos[1],x_distal->pos[2],pmj->pos[0],pmj->pos[1],pmj->pos[2]);
        if (dist < min_dist)
        {
            min_dist = dist;
            result = s;
        }
    }
    return result;
}

void Shocker::update_and_check_integrity ()
{
    // Check for null segments
    uint32_t ns = this->segment_list.size();
    for (uint32_t i = 0; i < ns; i++)
    {
        Segment *s = this->segment_list[i];
        s->calc_middle_point(s->middle_pos);
        if (s->length == 0.0)
        {
            fprintf(stderr,"[!] WARNING! Segment '%u' has a null length! src[%u]=(%g %g %g) || dest[%u]=(%g %g %g)\n",s->id,\
                    s->src->id,s->src->pos[0],s->src->pos[1],s->src->pos[2],\
                    s->dest->id,s->dest->pos[0],s->dest->pos[1],s->dest->pos[2]);
            exit(EXIT_FAILURE);
        }
    }
    // Update Min/Max LAT values
    Error *error = this->error;
    error->update_min_max_terminal_lat(this->segment_list);
}

double Shocker::update_minimum_distance (uint32_t &tosses, const double l_min)
{
    tosses++;
    if (tosses > NTOSS)  
    {
        //printf("Updating l_min -- %g to %g\n",l_min,l_min*FACTOR);
        tosses = 0;
        return l_min*FACTOR;
    }
    return l_min;
}

void Shocker::print_network_info ()
{
    //Graph *shocker_graph = convert_shocker_to_graph(this->node_list,this->segment_list);

    std::vector<double> segments;
    get_segment_length(this->segment_list,segments);
    
    double mean_segment_length, std_segment_length;
    calc_mean_std(segments,mean_segment_length,std_segment_length);
    write_vector_to_file(segments,this->params->output_dir + "/segments_length.dat");
        
    std::vector<double> angles;
    get_bifurcation_angles(this->segment_list,angles);
    
    double mean_biff_angle, std_biff_angle;
    calc_mean_std(angles,mean_biff_angle,std_biff_angle);
    write_vector_to_file(angles,this->params->output_dir + "/bifurcation_angle.dat");

    write_geometric_info_to_file(this->params->output_dir + "/network_info.txt",\
                    segments,mean_segment_length,std_segment_length,\
                    angles,mean_biff_angle,std_biff_angle);
    
    printf("------------------ NETWORK INFORMATION --------------------------------------------\n");
    printf("[INFO] Total number of segment = %lu\n",segments.size());
    printf("[INFO] Segment length = %.2lf +/- %.2lf mm\n",mean_segment_length,std_segment_length);
    printf("[INFO] Total number of bifurcations = %lu\n",angles.size());
    printf("[INFO] Bifurcation angle = %.2lf +/- %.2lf degrees\n",mean_biff_angle,std_biff_angle);
    printf("[INFO] Output files saved at: '%s'\n",this->params->output_dir.c_str());
    if (this->params->using_pmj_location) 
    {
        PMJ_Data *pmj_data = this->cloud->pmj_data;
        printf("[INFO] Number of PMJ's connected = %u/%lu [%g %%]\n",this->cloud->pmj_data->total_num_connected,this->cloud->pmj_data->pmjs.size(),(double)this->cloud->pmj_data->total_num_connected/(double)this->cloud->pmj_data->pmjs.size()*100.0);
        pmj_data->write_errors(this->params->output_dir + "/pmj_error.dat");
        pmj_data->write_errors_in_vtk(this->params->output_dir + "/pmj_error.vtk");

        this->error->calculate_electric_error(pmj_data);
        printf("[INFO] Max Error = %.2lf ms || Ref Min LAT = %.2lf ms || Ref Max LAT = %.2lf ms ||\n",this->error->max_lat_error,this->error->min_max_ref_lat[0],this->error->min_max_ref_lat[1]);                                                                                                           
        printf("                         || Aprox Min LAT = %.2lf ms || Aprox Max LAT = %.2lf ms ||\n",this->error->min_max_aprox_lat[0],this->error->min_max_aprox_lat[1]);
        printf("       RMSE = %.2lf ms || RRMSE = %.2lf %%\n",this->error->rmse,this->error->rrmse*100.0);    
        printf("       Epsilon 2ms = %.2lf %% || Epsilon = 5ms = %.2lf %%\n",this->error->epsilon_2ms*100.0,this->error->epsilon_5ms*100.0);       

        write_electric_info_to_file(this->params->output_dir + "/network_info.txt",\
                    this->error->max_lat_error,this->error->min_max_ref_lat[0],this->error->min_max_ref_lat[1],\
                    this->error->min_max_aprox_lat[0],this->error->min_max_aprox_lat[1],\
                    this->error->rmse,this->error->rrmse*100.0,\
                    this->error->epsilon_2ms*100.0,this->error->epsilon_5ms*100.0);
        printf("[INFO] Number of PMJs connected in the main loop: %u\n",this->error->counter_main_loop_connections);
        printf("[INFO] Number of PMJs connected after droping the LAT error tolerance: %u\n",this->error->counter_no_lat_error_tolerance_connections);
        printf("[INFO] Number of PMJs connected after droping the distance criterion: %u\n",this->error->counter_no_distance_criterion_connections);
        printf("[INFO] Number of PMJs connected after droping the distance criterion (geodesic): %u\n",this->error->counter_no_distance_criterion_connections_geodesic);
        printf("[INFO] Number of PMJs connected after droping the distance criterion (straight line): %u\n",this->error->counter_no_distance_criterion_connections_line);

    }
    printf("----------------------------------------------------------------------------------------\n");
    
    write_network_info_to_logfile(mean_segment_length,std_segment_length,segments.size(),\
                                  mean_biff_angle,std_biff_angle,angles.size());
    save_network_state();
}

void Shocker::load_network_state ()
{
    std::string filename = this->params->initial_network_name;
    std::cout << get_color_code(GRAY) << "[INFO] Making root using initial network '" << filename << "'" << std::endl;
    FILE *file = fopen(filename.c_str(),"r");
    if (!file)
    {
        fprintf(stderr,"[load_state] ERROR! Unable to open file '%s'!\n",filename.c_str());
        exit(EXIT_FAILURE);
    }
    uint32_t num_nodes;
    fscanf(file,"%u",&num_nodes);
    for (uint32_t i = 0; i < num_nodes; i++)
    {
        uint32_t id;
        double pos[3];
        double lat;
        int is_active;
        fscanf(file,"%u %lf %lf %lf %d",&id,&pos[0],&pos[1],&pos[2],&is_active);

        Node *n = new Node(id,pos,is_active);
        this->node_list.push_back(n);
    }
    uint32_t num_segments;
    fscanf(file,"%u",&num_segments);
    for (uint32_t i = 0; i < num_segments; i++)
    {
        uint32_t id;
        uint32_t edge[2];
        double lat;
        double diameter;
        bool mintree;
        fscanf(file,"%u %u %u %lf %lf %d",&id,&edge[0],&edge[1],&lat,&diameter,&mintree);
        Node *src = this->node_list[edge[0]];
        Node *dest = this->node_list[edge[1]];
        Segment *s = new Segment(id,diameter,lat,mintree,src,dest,nullptr,nullptr,nullptr);
        this->segment_list.push_back(s);
    }
    for (uint32_t i = 0; i < num_segments; i++)
    {
        int src_id, dest_id;
        fscanf(file,"%d %d",&src_id,&dest_id);
        this->segment_list[src_id]->parent = (dest_id != -1) ? this->segment_list[dest_id] : nullptr;
        fscanf(file,"%d %d",&src_id,&dest_id);
        this->segment_list[src_id]->right = (dest_id != -1) ? this->segment_list[dest_id] : nullptr;
        fscanf(file,"%d %d",&src_id,&dest_id);
        this->segment_list[src_id]->left = (dest_id != -1) ? this->segment_list[dest_id] : nullptr;
    }
    uint32_t num_terminals = 0;
    for (uint32_t i = 0; i < num_segments; i++)
        if (this->segment_list[i]->is_terminal()) num_terminals++;
    this->params->num_terminals = num_terminals;
    recalculate_length(this->segment_list);
    recalculate_local_activation_time(this->segment_list);
    fclose(file);
}

void Shocker::save_network_state ()
{
    std::cout << get_color_code(GRAY) << "[INFO] Saving network topology" << std::endl;
    std::string filename = this->params->output_dir + "/network_topology.txt";
    uint32_t num_points = this->node_list.size();
    uint32_t num_segments = this->segment_list.size();
    FILE *file = fopen(filename.c_str(),"w+");
    if (!file)
    {
        fprintf(stderr,"[shocker] ERROR! Opening file '%s'\n",filename.c_str());
        exit(EXIT_FAILURE);
    }
    fprintf(file,"%u\n",num_points);
    for (uint32_t i = 0; i < num_points; i++)
        fprintf(file,"%u %g %g %g %d\n",this->node_list[i]->id,\
                                this->node_list[i]->pos[0],this->node_list[i]->pos[1],this->node_list[i]->pos[2],\
                                static_cast<int>(this->node_list[i]->is_active));
    fprintf(file,"%u\n",num_segments);
    for (uint32_t i = 0; i < num_segments; i++)
    {
        uint32_t cur_seg_id = this->segment_list[i]->id;
        uint32_t src_id = this->segment_list[i]->src->id;
        uint32_t dest_id = this->segment_list[i]->dest->id;
        double lat = this->segment_list[i]->lat;
        double diameter = this->segment_list[i]->diameter;
        bool is_mintree = this->segment_list[i]->mintree;
        fprintf(file,"%u %u %u %g %g %d\n",cur_seg_id,src_id,dest_id,lat,diameter,static_cast<int>(is_mintree));
    }
    for (uint32_t i = 0; i < num_segments; i++)
    {
        int cur_id = this->segment_list[i]->id;
        int parent_id = (this->segment_list[i]->parent) ? this->segment_list[i]->parent->id : -1;
        int left_id = (this->segment_list[i]->left) ? this->segment_list[i]->left->id : -1;
        int right_id = (this->segment_list[i]->right) ? this->segment_list[i]->right->id : -1;
        fprintf(file,"%d %d\n",cur_id,parent_id);
        fprintf(file,"%d %d\n",cur_id,left_id);
        fprintf(file,"%d %d\n",cur_id,right_id);
    }
    fclose(file);
}

void Shocker::reset_network_state ()
{
    for (uint32_t i = 0; i < this->segment_list.size(); i++) 
        delete this->segment_list[i];
    this->segment_list.clear();
    for (uint32_t i = 0; i < this->node_list.size(); i++) 
        delete this->node_list[i];
    this->node_list.clear();
}

void Shocker::build_minimum_network_maps (std::map<uint32_t,uint32_t> &nodes_map, std::map<uint32_t,uint32_t> &segments_map)
{
    std::map<uint32_t,uint32_t>::iterator it;
    uint32_t nn = this->node_list.size();
    uint32_t ns = this->segment_list.size();
    uint32_t counter_points = 0;
    uint32_t counter_segments = 0;
    for (uint32_t i = 0; i < ns; i++)
    {
        Segment *s = this->segment_list[i];

        // We only check the segments which are tagged as part of the minimum network
        if (s->is_mintree())
        {
            Node *src = s->src;
            Node *dest = s->dest;

            it = nodes_map.find(src->id);
            if (it == nodes_map.end())     // Point not inserted yet ...
            {
                nodes_map.insert( std::make_pair(src->id,counter_points) );
                counter_points++;    
            }
            it = nodes_map.find(dest->id);
            if (it == nodes_map.end())     // Point not inserted yet ...
            {
                nodes_map.insert( std::make_pair(dest->id,counter_points) );
                counter_points++;    
            }
            it = segments_map.find(s->id);
            if (it == segments_map.end())   // Segment not inserted yet ...
            {
                segments_map.insert( std::make_pair(s->id,counter_segments) );
                counter_segments++;
            }
        }
    }
}

void Shocker::build_minimum_network_topology (std::map<uint32_t,uint32_t> nodes_map, std::map<uint32_t,uint32_t> segments_map,\
                                std::vector<Node> &min_network_nodes, std::vector<Segment> &min_network_segments)
{
    std::map<uint32_t,uint32_t>::iterator it;
    // Initialize the network points
    min_network_nodes.assign(nodes_map.size(),Node());
    for (it = nodes_map.begin(); it != nodes_map.end(); ++it)
    {
        uint32_t full_network_id = it->first;
        uint32_t min_network_id = it->second;
        Node *full_network_node = this->node_list[full_network_id]; 
        Node *min_network_node = &min_network_nodes[min_network_id];

        min_network_node->setNode(full_network_node);
        min_network_node->setId(min_network_id);
    }
    // Initialize the network segments
    min_network_segments.assign(segments_map.size(),Segment());
    for (it = segments_map.begin(); it != segments_map.end(); ++it)
    {
        uint32_t full_network_id = it->first;
        uint32_t min_network_id = it->second;

        Node *full_network_src = this->segment_list[full_network_id]->src;
        Node *full_network_dest = this->segment_list[full_network_id]->dest;
        
        uint32_t min_network_src_id = nodes_map[full_network_src->id];
        uint32_t min_network_dest_id = nodes_map[full_network_dest->id];

        Node *min_network_src = &min_network_nodes[min_network_src_id];
        Node *min_network_dest = &min_network_nodes[min_network_dest_id];
        
        Segment *min_network_segment = &min_network_segments[min_network_id];
        min_network_segment->id = min_network_id;
        min_network_segment->src = min_network_src;
        min_network_segment->dest = min_network_dest;
        min_network_segment->diameter = this->segment_list[full_network_id]->diameter;
        min_network_segment->length = this->segment_list[full_network_id]->length;
        min_network_segment->lat = this->segment_list[full_network_id]->lat;
        min_network_segment->mintree = this->segment_list[full_network_id]->mintree;        
    }
    // Set the segment pointers
    for (it = segments_map.begin(); it != segments_map.end(); ++it)
    {
        uint32_t full_network_id = it->first;
        uint32_t min_network_id = it->second;

        Segment *full_network_parent = this->segment_list[full_network_id]->parent;
        Segment *full_network_left = this->segment_list[full_network_id]->left;
        Segment *full_network_right = this->segment_list[full_network_id]->right;
        
        int min_network_parent_id = (full_network_parent && full_network_parent->is_mintree()) ? segments_map[full_network_parent->id] : -1;
        int min_network_left_id = (full_network_left && full_network_left->is_mintree()) ? segments_map[full_network_left->id] : -1;
        int min_network_right_id = (full_network_right && full_network_right->is_mintree()) ? segments_map[full_network_right->id] : -1;
        
        Segment *min_network_parent = (min_network_parent_id != -1) ? &min_network_segments[min_network_parent_id] : nullptr;
        Segment *min_network_left = (min_network_left_id != -1) ? &min_network_segments[min_network_left_id] : nullptr;
        Segment *min_network_right = (min_network_right_id != -1) ? &min_network_segments[min_network_right_id] : nullptr;

        min_network_segments[min_network_id].parent = min_network_parent;
        min_network_segments[min_network_id].left = min_network_left;
        min_network_segments[min_network_id].right = min_network_right;
    }
}

void Shocker::recalculate_errors ()
{
    PMJ_Data *pmj_data = this->cloud->pmj_data;
    double his_offset = this->params->his_offset;
    for (uint32_t i = 0; i < pmj_data->pmjs.size(); i++)
    {
        PMJPoint *pmj = pmj_data->pmjs[i];
        Segment *iterm = search_active_terminal_segment(pmj);
        pmj->aprox_value = iterm->calc_terminal_local_activation_time() + his_offset; 
        pmj->error = (pmj->ref_value - pmj->aprox_value);
    }
}

bool Shocker::attempt_diameter_adjustment ()
{
    PMJ_Data *pmj_data = this->cloud->pmj_data;
    for (uint32_t i = 0; i < pmj_data->pmjs.size(); i++)
    {
        PMJPoint *pmj = pmj_data->pmjs[i];
        Segment *iterm = search_active_terminal_segment(pmj);
        Segment *tmp = iterm;
        double dist = 0.0;
        while (!tmp->is_bifurcation())
        {
            dist += tmp->length;
            tmp = tmp->parent;
        }
        double t1 = tmp->lat;
        double t2 = pmj->ref_value;
        double max_diam = tmp->diameter;
        double cv_target = (dist / (t2-t1))*UM_MS_TO_M_S;
        double diam_target = tmp->calc_diameter(cv_target);
        // Catch an error when the terminal branch has a larget diameter larger than its parent
        if (diam_target > max_diam)
        {
            std::cout << get_color_code(RED) << "[-] ERROR! It is impossible to adjust the diameter for PMJ " << pmj->id << "! Non-decreasing diameter error!" << std::endl;
            exit(EXIT_FAILURE);    
        }
        std::cout << get_color_code(YELLOW) << "[PMJ " << pmj->id << "] Updated diameter from " << iterm->diameter << "um to " << diam_target << "um" << std::endl;
        std::cout << get_color_code(YELLOW) << "         Updated CV from " << iterm->calc_propagation_velocity() << "m/s to " << cv_target << "m/s" << std::endl; 
        tmp = iterm;
        while (!tmp->is_bifurcation())
        {
            tmp->diameter = diam_target;
            tmp = tmp->parent;
        }
    }
    recalculate_local_activation_time(this->segment_list);
    recalculate_errors();
    std::cout << get_color_code(GREEN);
    pmj_data->print();
    return true;
}

std::vector<std::vector<double>> Shocker::get_nodes_coordinates ()
{
    std::vector<std::vector<double>> result;
    for (uint32_t i = 0; i < this->node_list.size(); i++)
    {
        Node *node = this->node_list[i];
        result.push_back({node->pos[0],node->pos[1],node->pos[2]});
    }
    return result;
}

void Shocker::write_full_network_to_vtk ()
{
    Logger *logger = this->logger;
    logger->write_full_network_to_vtk(this->node_list,this->segment_list,\
                                    this->params->output_dir,this->params->his_offset,this->params->num_terminals);
    logger->write_full_network_monoalg_to_vtk(this->node_list,this->segment_list,\
                                    this->params->output_dir,this->params->his_offset,this->params->num_terminals);
}

void Shocker::write_minimum_network_to_vtk ()
{
    Logger *logger = this->logger;
    logger->write_minimum_network_to_vtk(this->node_list,this->segment_list,\
                                    this->params->output_dir,this->params->his_offset,this->params->num_terminals);
    logger->write_minimum_network_monoalg_to_vtk(this->node_list,this->segment_list,\
                                    this->params->output_dir,this->params->his_offset,this->params->num_terminals);
}

Segment* Shocker::build_branch (Segment *iconn, const double pos[], double &branch_size)
{
    uint32_t nn = this->node_list.size();
    uint32_t ns = this->segment_list.size();
    double d = iconn->diameter;
    branch_size = 0.0;

    start_stop_watch(&geodesic_pathway_time);

    Surface_Data *surface_data = this->cloud->surface_data;
    std::vector<CloudPoint> geodesic_points = surface_data->compute_geodesic_pathway(iconn->dest->pos,pos,false,iconn->middle_pos);
    
    total_time_geodesic_pathway += stop_stop_watch(&geodesic_pathway_time);

    bool has_duplicates = check_duplicates(this->node_list,geodesic_points);
    bool has_more_points = (geodesic_points.size() > 1);
    if (has_duplicates || !has_more_points)
    {
        //if (has_duplicates)
        //    printf("[build_branch] Duplicate points!\n");
        //if (!has_more_points)
        //    printf("[build_branch] Invalid geodesic pathway!\n");
        return nullptr;
    }
        
    //write_cloud_array_to_file("outputs/geodesic_path.vtk",geodesic_points);

    // Insert the new nodes
    for (uint32_t i = 0; i < geodesic_points.size(); i++)
    {
        double *node_pos = geodesic_points[i].pos;
        Node *node = new Node(nn+i,node_pos);
        this->node_list.push_back(node);
    }
    Node *biff_node = this->node_list[nn];

    // Create 'ibiff'
    Segment *ibiff = new Segment(ns,d,0,false,iconn->src,biff_node,nullptr,nullptr,iconn->parent);
    this->segment_list.push_back(ibiff);
    
    // Adjust 'iconn' pointers
    iconn->parent = ibiff;
    iconn->src = biff_node;
    ibiff->right = iconn;

    // Create 'inew' branch
    for (uint32_t i = 0; i < geodesic_points.size()-1; i++)
    {
        Segment *segment = new Segment(ns+i+1,d,0,false,this->node_list[nn+i],this->node_list[nn+i+1],nullptr,nullptr,nullptr);
        this->segment_list.push_back(segment);  
        branch_size += segment->length;
    }
    
    // Adjust 'inew' parent pointer
    Segment *inew = this->segment_list[ns+1];
    inew->parent = ibiff;

    // Adjust 'ibiff' pointers
    ibiff->left = inew;

    // Adjust 'inew' branch pointers
    for (uint32_t i = 0; i < geodesic_points.size()-2; i++)
    {
        inew->right = this->segment_list[ns+i+2];
        inew->right->parent = inew;
        inew = inew->right;
    }

    // Update the network parameters
    this->params->num_terminals++;
    recalculate_length(this->segment_list);
    recalculate_local_activation_time(this->segment_list);
    //recalculate_length_branch(inew,ibiff);
    //recalculate_local_activation_time_branch(inew,ibiff);

    //update_and_check_integrity();

    return inew;
}

Segment* Shocker::prune_branch (Segment *inew)
{
    Segment *aux = nullptr;
    Segment *tmp = inew;
    uint32_t counter_segments = 0;

    // Unset the pointers of the target branch
    while (!tmp->is_bifurcation())
    {
        aux = tmp->parent;

        tmp->parent = nullptr;
        tmp->left = nullptr;
        tmp->right = nullptr;
        tmp->src = nullptr;
        tmp->dest = nullptr;

        tmp = aux;
        counter_segments++;
    }

    // Update the pointers of the 'iconn'
    tmp->left = nullptr;
    tmp->right->src = tmp->src;
    tmp->right->parent = tmp->parent;

    // Remove the segment of the branch
    while (counter_segments > 0)
    {
        eliminate_node_from_list(this->node_list,this->node_list.back());
        eliminate_segment_from_list(this->segment_list,this->segment_list.back());
        counter_segments--;
    }
    eliminate_node_from_list(this->node_list,this->node_list.back());
    eliminate_segment_from_list(this->segment_list,this->segment_list.back());

    // Update the network parameters
    this->params->num_terminals--;
    recalculate_length(this->segment_list);
    recalculate_local_activation_time(this->segment_list);

    return tmp;
}

bool Shocker::check_duplicates (std::vector<Node*> n_list, std::vector<CloudPoint> geo_points)
{
    start_stop_watch(&check_duplicates_time);

    for (uint32_t i = 0; i < geo_points.size(); i++)
    {
        CloudPoint cp = geo_points[i];
        for (uint32_t j = 0; j < n_list.size(); j++)
        {
            Node *n = this->node_list[j];
            double dist = euclidean_norm(n->pos[0],n->pos[1],n->pos[2],cp.pos[0],cp.pos[1],cp.pos[2]);
            if (dist < 1.0e-08)
                return true;
        }
    }
    
    total_time_check_duplicates += stop_stop_watch(&check_duplicates_time);

    return false;
}

void Shocker::extract_connected_pmjs ()
{
    std::vector<PMJPoint*> pmjs_connected;
    PMJ_Data *pmj_data = this->cloud->pmj_data;
    for (uint32_t i = 0; i < pmj_data->pmjs.size(); i++) 
    {
        PMJPoint *pmj = pmj_data->pmjs[i];
        if (pmj->connected)
        {
            pmjs_connected.push_back(pmj);
        }
    }

    for (uint32_t i = 0; i < pmjs_connected.size(); i++)
    {
        PMJPoint *pmj = pmjs_connected[i];
        printf("%g %g %g %g\n",pmj->pos[0],pmj->pos[1],pmj->pos[2],pmj->ref_value);
    }

    //FILE *file = fopen("outputs/pmjs.vtk","w+");
    //fprintf(file,"# vtk DataFile Version 4.1\n");
    //fprintf(file,"vtk output\n");
    //fprintf(file,"ASCII\n");
    //fprintf(file,"DATASET POLYDATA\n");
    //fprintf(file,"POINTS %u float\n",pmjs_connected.size());
    //for (uint32_t i = 0; i < pmjs_connected.size(); i++)
    //    fprintf(file,"%g %g %g\n",pmjs_connected[i]->pos[0],pmjs_connected[i]->pos[1],pmjs_connected[i]->pos[2]);
    //fprintf(file,"VERTICES %u %u\n",pmjs_connected.size(),pmjs_connected.size()*2);
    //for (uint32_t i = 0; i < pmjs_connected.size(); i++)
    //    fprintf(file,"1 %u\n",i);
    //fclose(file);
}

void Shocker::write_network_info_to_logfile (const double mean_segment_length, const double std_segment_length, const uint32_t ns,\
                                             const double mean_biff_angle, const double std_biff_angle, const uint32_t nb)
{
    FILE *log_file = this->params->log_file;
    fprintf(log_file,"------------------ NETWORK INFORMATION --------------------------------------------\n");
    fprintf(log_file,"[INFO] Total number of segment = %lu\n",ns);
    fprintf(log_file,"[INFO] Segment length = %.2lf +/- %.2lf mm\n",mean_segment_length,std_segment_length);
    fprintf(log_file,"[INFO] Total number of bifurcations = %lu\n",nb);
    fprintf(log_file,"[INFO] Bifurcation angle = %.2lf +/- %.2lf degrees\n",mean_biff_angle,std_biff_angle);
    fprintf(log_file,"[INFO] Output files saved at: '%s'\n",this->params->output_dir.c_str());
    if (this->params->using_pmj_location) 
    {
        PMJ_Data *pmj_data = this->cloud->pmj_data;
        fprintf(log_file,"[INFO] Number of PMJ's connected = %u/%lu [%g %%]\n",this->cloud->pmj_data->total_num_connected,\
                                                                            this->cloud->pmj_data->pmjs.size(),\
                                                                            (double)this->cloud->pmj_data->total_num_connected/(double)this->cloud->pmj_data->pmjs.size()*100.0);

        fprintf(log_file,"[INFO] Max Error = %.2lf ms || Ref Min LAT = %.2lf ms || Ref Max LAT = %.2lf ms ||\n",this->error->max_lat_error,this->error->min_max_ref_lat[0],this->error->min_max_ref_lat[1]);                                                                                                           
        fprintf(log_file,"                         || Aprox Min LAT = %.2lf ms || Aprox Max LAT = %.2lf ms ||\n",this->error->min_max_aprox_lat[0],this->error->min_max_aprox_lat[1]);
        fprintf(log_file,"       RMSE = %.2lf ms || RRMSE = %.2lf %%\n",this->error->rmse,this->error->rrmse*100.0);    
        fprintf(log_file,"       Epsilon 2ms = %.2lf %% || Epsilon = 5ms = %.2lf %%\n",this->error->epsilon_2ms*100.0,this->error->epsilon_5ms*100.0);       

        fprintf(log_file,"[INFO] Number of PMJs connected in the main loop: %u\n",this->error->counter_main_loop_connections);
        fprintf(log_file,"[INFO] Number of PMJs connected after droping the LAT error tolerance: %u\n",this->error->counter_no_lat_error_tolerance_connections);  
        fprintf(log_file,"[INFO] Number of PMJs connected after droping the distance criterion: %u\n",this->error->counter_no_distance_criterion_connections);  
        fprintf(log_file,"[INFO] Number of PMJs connected after droping the distance criterion (geodesic): %u\n",this->error->counter_no_distance_criterion_connections_geodesic);
        fprintf(log_file,"[INFO] Number of PMJs connected after droping the distance criterion (straight line): %u\n",this->error->counter_no_distance_criterion_connections_line);
    }
    fprintf(log_file,"----------------------------------------------------------------------------------------\n");
}

bool Shocker::refinement_attempt_diameter_adjustment ()
{
    std::cout << get_color_code(YELLOW) << "[shocker] Refinement of the solution by adjusting the terminal branches diameter" << std::endl;

    PMJ_Data *pmj_data = this->cloud->pmj_data;
    double his_offset = this->params->his_offset;
    Segment *iterm = nullptr;
    for (uint32_t i = 0; i < pmj_data->pmjs.size(); i++)
    {
        PMJPoint *pmj = pmj_data->pmjs[i];
        iterm = search_active_terminal_segment(pmj);
        
        double branch_size = iterm->length;
        Segment *tmp = iterm;
        while (!tmp->is_bifurcation())
        {
            branch_size += tmp->length;
            tmp = tmp->parent;
        }

        double t1 = tmp->calc_terminal_local_activation_time() + his_offset;
        double t2 = pmj->ref_value;
        if (pmj->error > 0)
        {
            std::cout << get_color_code(CYAN) << "[DEBUG] We can adjust PMJ point " << pmj->id << " has PMJ error equals to " << pmj->error << std::endl;
            std::cout << get_color_code(CYAN) << "[DEBUG] ibiff LAT is " << t1 << " and iterm LAT is " << t2 << std::endl;
            
            double cv_original = tmp->calc_propagation_velocity();
            double cv_target = (branch_size / (t2-t1))*UM_MS_TO_M_S;
            double d_target = tmp->calc_diameter(cv_target);
            double d_original = iterm->diameter;
            if (cv_target >= 1.0 && cv_target <= 4.0 && d_target <= d_original)
            {
                std::cout << get_color_code(GREEN) << "[shocker] We can adjust the branch diameter for PMJ point " << pmj->id << std::endl;
                std::cout << get_color_code(GREEN) << "\t[shocker] Decreasing diameter from " << d_original << " to " << d_target << std::endl;
                std::cout << get_color_code(GREEN) << "\t[shocker] Decreasing CV from " << cv_original << " to " << cv_target << std::endl;
                
                tmp = iterm;
                while (!tmp->is_bifurcation())
                {
                    tmp->diameter = d_target;
                    tmp = tmp->parent;
                }
            }
            else
            {
                std::cout << get_color_code(RED) << "[shocker] ERROR! The adjustable diameter for PMJ point " << pmj->id << " will be unphysiological!" << std::endl;
                std::cout << get_color_code(RED) << "\t[shocker] Decreasing diameter from " << d_original << " to " << d_target << std::endl;
                std::cout << get_color_code(RED) << "\t[shocker] Decreasing CV from " << cv_original << " to " << cv_target << std::endl;
            }
        }
    }
    recalculate_length(this->segment_list);
    recalculate_local_activation_time(this->segment_list);
    recalculate_errors();

    write_minimum_network_to_vtk();

    return true;
}