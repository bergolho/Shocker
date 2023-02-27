#include "cloud.h"

Cloud::Cloud (CloudConfig *cloud_config, PMJConfig *pmj_config, std::string output_dir, const double root_pos[])
{
    this->cloud_data = new Cloud_Data(cloud_config,output_dir,root_pos);
    this->pmj_data = new PMJ_Data(pmj_config,output_dir);
    this->surface_data = new Surface_Data(cloud_config->surface_filename);
    this->write(output_dir);
}

Cloud::~Cloud ()
{
    if (this->cloud_data)
        delete this->cloud_data;
    if (this->pmj_data)
        delete this->pmj_data;
    if (this->surface_data)
        delete this->surface_data;
}

void Cloud::print ()
{
    this->cloud_data->print();
    if (this->pmj_data->has_pmjs())
        this->pmj_data->print();
}

void Cloud::write (std::string output_dir)
{
    this->cloud_data->write(output_dir);
    if (this->pmj_data->has_pmjs())
        this->pmj_data->write(output_dir);
}

void Cloud::tag_pmjs_to_cloud ()
{
    uint32_t num_cloud_points = this->cloud_data->points.size();
    uint32_t num_pmjs = this->pmj_data->pmjs.size();
    for (uint32_t i = 0; i < num_cloud_points; i++)
    {
        CloudPoint *c_point = this->cloud_data->points[i];
        double closest_dist = __DBL_MAX__;
        PMJPoint *closest = nullptr;
    
        for (uint32_t j = 0; j < num_pmjs; j++)
        {
            PMJPoint *p_point = this->pmj_data->pmjs[j];
            double dist = euclidean_norm(c_point->pos[0],c_point->pos[1],c_point->pos[2],\
                                        p_point->pos[0],p_point->pos[1],p_point->pos[2]);
            if (dist < closest_dist)
            {
                closest_dist = dist;
                closest = p_point;
            }
        }
        c_point->closest_pmj = closest;
    }
}

Cloud_Data::Cloud_Data (CloudConfig *cloud_config, std::string output_dir, const double root_pos[])
{
    this->cloud_filename = cloud_config->cloud_filename;
    this->phi = cloud_config->phi;
    this->rand_offset = cloud_config->rand_offset;
    this->cur_id = 0;
    
    bool sucess = read_cloud_from_vtk();
    if (sucess) 
    {
        printf("[cloud] Cloud of points was sucessfully loaded from: '%s'!\n",cloud_config->cloud_filename.c_str());
    }
    else
    {
        fprintf(stderr,"[-] ERROR! Cannot read cloud of points file: '%s'!\n",cloud_config->cloud_filename.c_str());
        exit(EXIT_FAILURE);
    }
    
    // PRE-PROCESSING
    pre_process(root_pos);
}

Cloud_Data::~Cloud_Data ()
{
    for (uint32_t i = 0; i < this->points.size(); i++)
        delete this->points[i];
    this->points.clear();
}

CloudPoint* Cloud_Data::sort_point ()
{
    uint32_t prev_id = this->cur_id;
    uint32_t offset = rand() % this->rand_offset + 1;
    this->cur_id = (this->cur_id + offset) % (this->points.size());
    while (this->points[this->cur_id]->taken)
    {
        offset = rand() % this->rand_offset;
        this->cur_id = (this->cur_id + offset) % (this->points.size());
    }
    
    if (this->cur_id < prev_id)
    {
        this->counter_cloud_passes++;
        //printf("Pass through the whole cloud of points once!\n");
    }
    //printf("cur_id: %u/%u\n",this->cur_id,this->points.size());
    return this->points[this->cur_id];
}

void Cloud_Data::print ()
{
    printf("[cloud] Cloud points filename: '%s'\n",this->cloud_filename.c_str());
    printf("[cloud] rand_offset = %g\n",this->rand_offset);
    printf("[cloud] phi = %g\n",this->phi);
    printf("[cloud] Number of points: %lu\n",this->points.size());
    for (uint32_t i = 0; i < this->points.size(); i++)
        this->points[i]->print();
}

void Cloud_Data::write (std::string output_dir)
{
    assert(this->points.size() > 0);

    uint32_t np = this->points.size();
    std::string name = output_dir + "/cloud_points_um.vtk";
    FILE *file = fopen(name.c_str(),"w+");

    fprintf(file,"# vtk DataFile Version 4.1\n");
    fprintf(file,"vtk output\n");
    fprintf(file,"ASCII\n");
    fprintf(file,"DATASET POLYDATA\n");
    fprintf(file,"POINTS %u float\n",np);
    for (uint32_t i = 0; i < np; i++)
        fprintf(file,"%g %g %g\n",this->points[i]->pos[0],this->points[i]->pos[1],this->points[i]->pos[2]);
        
    fprintf(file,"VERTICES %u %u\n",np,np*2);
    for (uint32_t i = 0; i < np; i++)
        fprintf(file,"1 %u\n",i);

    fprintf(file,"POINT_DATA %u\n",np);
    fprintf(file,"FIELD FieldData 1\n");
    fprintf(file,"ClosestPMJ 1 %u float\n",np);
    for (uint32_t i = 0; i < np; i++)
    {
        if (this->points[i]->closest_pmj)
            fprintf(file,"%u\n",this->points[i]->closest_pmj->id);
        else
            fprintf(file,"-1\n");
    }

    fclose(file);
}

void Cloud_Data::pre_process (const double root_pos[])
{
    remap(root_pos);
    filter();
}

void Cloud_Data::remap (const double root_pos[])
{
    std::vector<CloudPoint*> remapped_points;
    uint32_t np = this->points.size();
    //char filename[200];

    // Remapping the points
    // Use a growing sphere centered at the root position to get the nearest points
    uint32_t counter = 0;
    double region_radius = MIN_REGION_RADIUS;
    while (remapped_points.size() < np)
    {
        for (uint32_t i = 0; i < this->points.size(); i++)
        {
            CloudPoint *p = this->points[i];
            double dist = euclidean_norm(p->pos[0],p->pos[1],p->pos[2],root_pos[0],root_pos[1],root_pos[2]);
            if (dist < region_radius)
            {
                remapped_points.push_back(p);
                this->points.erase(this->points.begin()+i); // Delete the unmapped point
            }    
        }
        region_radius += REGION_RADIUS_OFFSET;

        //printf("Counter = %u --> radius = %g\n",counter,region_radius);
        //sprintf(filename,"cloud_figure/remapped_cloud_%u.vtk",counter);
        //write_cloud_array_to_file(filename,remapped_points);
        counter++;
    }

    for (uint32_t i = 0; i < remapped_points.size(); i++)
    {
        CloudPoint *p = remapped_points[i];
        this->points.push_back(p);
    }
    remapped_points.clear();
}

void Cloud_Data::filter ()
{
    std::vector<CloudPoint*> filtered_points;

    // Number of points
    const uint32_t np = this->points.size();
    const uint32_t rnp = np*this->phi;

    // Use a bitmask to set the points to be taken
    uint32_t counter = 0;
    std::vector<bool> bitmask;
    bitmask.assign(np,false);
    while (counter < rnp)
    {
        int k = rand() % bitmask.size();
        if (!bitmask[k])
        {
            bitmask[k] = true;
            counter++;
            if (counter == rnp) break;
        }
    }

    // Filter the points using the bitmask
    for (uint32_t i = 0; i < np; i++)
    {
        CloudPoint *p = this->points[i];
        if (bitmask[i]) 
            filtered_points.push_back(p);
        else
            delete this->points[i];
    }
    this->points.clear();
    
    // Fill the points with the filtered ones
    for (uint32_t i = 0; i < filtered_points.size(); i++)
    {
        CloudPoint *p = filtered_points[i];
        this->points.push_back(p);
    }

    //char filename[200];
    //sprintf(filename,"cloud_figure/filtered_cloud.vtk");
    //write_cloud_array_to_file(filename,this->points); 

    //exit(EXIT_SUCCESS);   
}

bool Cloud_Data::read_cloud_from_vtk ()
{
    char str[500];
    double pos[3];
    uint32_t vertex[2];
    uint32_t num_points;

    FILE *file = fopen(this->cloud_filename.c_str(),"r");
    if (!file) return false;

    while (fscanf(file,"%s",str) != EOF)
        if (strcmp(str,"POINTS") == 0) break;
    
    fscanf(file,"%u %s",&num_points,str);
    for (uint32_t i = 0; i < num_points; i++)
    {    
        fscanf(file,"%lf %lf %lf",&pos[0],&pos[1],&pos[2]);
        CloudPoint *p = new CloudPoint(i,pos,0,false);
        points.push_back(p);
    }

    while (fscanf(file,"%s",str) != EOF)
        if (strcmp(str,"VERTICES") == 0) break;

    int trash[2];
    fscanf(file,"%d %d",&trash[0],&trash[1]);

    for (uint32_t i = 0; i < num_points; i++)
    {
        fscanf(file,"%u %u",&vertex[0],&vertex[1]);
    }

    // Check if there is POINT_DATA
    int is_eof = fscanf(file,"%s",str);
    if ( is_eof != EOF && strcmp(str,"POINT_DATA") == 0 )
    {
        while (fscanf(file,"%s",str) != EOF)
            if (strcmp(str,"float") == 0) break;
        
        // Read the LAT section
        for (uint32_t i = 0; i < num_points; i++)
        {
            double value;
            fscanf(file,"%lf",&value);
            points[i]->setLAT(value);
        }
    }
    fclose(file);

    return true;
}

void Cloud_Data::calc_distal_root_position (const double x_prox[], const double l_d, double x_dist[])
{
    uint32_t index = 0;
    uint32_t counter = 0;
    uint32_t iterations = 0;
    bool is_root_ok = false;
    double l_min = sqrt( l_d * l_d );   // Minimum length
    CloudPoint *p = nullptr;
    // Sort a point until we surpass the minimum length
    while (!is_root_ok)
    {
        p = sort_point();
        double l = euclidean_norm(x_prox[0],x_prox[1],x_prox[2],p->pos[0],p->pos[1],p->pos[2]);

        if (l >= l_min)
            is_root_ok = true;
        else
            counter++;
            
        if (counter > 8) l_min *= 0.95;

        iterations++;   
    }
    memcpy(x_dist,p->pos,sizeof(double)*3);
    printf("[cloud] Root segment was set in %u iterations.\n", iterations);
    printf("[cloud] Root distal position set to: (%g %g %g)\n",p->pos[0],p->pos[1],p->pos[2]);
}

PMJ_Data::PMJ_Data (PMJConfig *config, std::string output_dir)
{
    if (config->using_pmj)
    {
        this->total_num_connected = 0;
        this->connection_rate = config->connection_rate;
        this->lat_error_tolerance = config->lat_error_tolerance;
        this->pmj_filename = config->location_filename;

        bool sucess = read_pmjs_from_vtk();
        if (sucess) 
        {
            printf("[pmj] PMJ location file was sucessfully loaded from: '%s'!\n",this->pmj_filename.c_str());
        }
        else
        {
            fprintf(stderr,"[pmj] ERROR! Cannot read PMJ location file: '%s'!\n",config->location_filename.c_str());
            exit(EXIT_FAILURE);
        }
        // Sort the PMJs by their LAT
        std::sort(this->pmjs.begin(),this->pmjs.end(),compareByLAT);
    }
    else
    {
        this->total_num_connected = 0;
        this->connection_rate = __UINT32_MAX__;
        this->lat_error_tolerance = 2.0;
        this->pmj_filename = "";
    }
    //print();
}

PMJ_Data::~PMJ_Data ()
{
    for (uint32_t i = 0; i < this->pmjs.size(); i++)
        delete this->pmjs[i];
}

void PMJ_Data::reset ()
{
    this->total_num_connected = 0;
    for (uint32_t i = 0; i < this->pmjs.size(); i++)
    {
        PMJPoint *pmj = this->pmjs[i];
        pmj->connected = false;
        pmj->aprox_value = 0.0;
        pmj->error = 0.0;
    }
}

bool PMJ_Data::has_pmjs ()
{
    return (this->pmjs.size() > 0);
}

void PMJ_Data::print ()
{
    assert(this->pmjs.size() > 0);

    printf("[pmj] Total number of connected PMJs = %u/%u\n",this->total_num_connected,this->pmjs.size());
    printf("[pmj] Connection rate = %u\n",this->connection_rate);
    printf("[pmj] LAT error tolerance = %g\n",this->lat_error_tolerance);
    printf("[pmj] PMJ location filename = %s\n",this->pmj_filename.c_str());
    for (uint32_t i = 0; i < this->pmjs.size(); i++)
        this->pmjs[i]->print();
}

bool PMJ_Data::read_pmjs_from_vtk ()
{
    char str[500];
    double pos[3];
    uint32_t vertex[2];
    uint32_t num_points;

    FILE *file = fopen(this->pmj_filename.c_str(),"r");
    if (!file) return false;

    while (fscanf(file,"%s",str) != EOF)
        if (strcmp(str,"POINTS") == 0) break;
    
    fscanf(file,"%u %s",&num_points,str);
    for (uint32_t i = 0; i < num_points; i++)
    {
        fscanf(file,"%lf %lf %lf",&pos[0],&pos[1],&pos[2]);
        PMJPoint *pmj = new PMJPoint(i,pos);
        this->pmjs.push_back(pmj);
    }

    while (fscanf(file,"%s",str) != EOF)
        if (strcmp(str,"VERTICES") == 0) break;

    int trash[2];
    fscanf(file,"%d %d",&trash[0],&trash[1]);

    for (uint32_t i = 0; i < num_points; i++)
    {
        uint32_t vertex[2];
        fscanf(file,"%u %u",&vertex[0],&vertex[1]);
    }

    // Check if there is POINT_DATA section in the file
    int is_eof = fscanf(file,"%s",str);
    if ( is_eof != EOF && strcmp(str,"POINT_DATA") == 0 )
    {
        while (fscanf(file,"%s",str) != EOF)
            if (strcmp(str,"float") == 0) break;
        
        // Read the reference LAT section
        for (uint32_t i = 0; i < num_points; i++)
        {
            double value;
            fscanf(file,"%lf",&value);
            pmjs[i]->ref_value = value;
        }
    }  
    fclose(file);

    return true;
}

void PMJ_Data::write (std::string output_dir)
{
    assert(this->pmjs.size() > 0);

    uint32_t np = this->pmjs.size();
    std::string name = output_dir + "/pmj_points_um.vtk";
    FILE *file = fopen(name.c_str(),"w+");

    fprintf(file,"# vtk DataFile Version 4.1\n");
    fprintf(file,"vtk output\n");
    fprintf(file,"ASCII\n");
    fprintf(file,"DATASET POLYDATA\n");
    fprintf(file,"POINTS %u float\n",np);
    for (uint32_t i = 0; i < np; i++)
        fprintf(file,"%g %g %g\n",this->pmjs[i]->pos[0],this->pmjs[i]->pos[1],this->pmjs[i]->pos[2]);
        
    fprintf(file,"VERTICES %u %u\n",np,np*2);
    for (uint32_t i = 0; i < np; i++)
        fprintf(file,"1 %u\n",i);
    
    fprintf(file,"POINT_DATA %u\n",np);
    fprintf(file,"FIELD FieldData 1\n");
    fprintf(file,"Ref_LAT 1 %u float\n",np);
    for (uint32_t i = 0; i < np; i++)
        fprintf(file,"%g\n",this->pmjs[i]->ref_value);
    
    fclose(file);
}

void PMJ_Data::write_errors (std::string filename)
{
    FILE *file = fopen(filename.c_str(),"w+");
    fprintf(file,"[INFO] Number of connected PMJs: %u/%u\n",this->total_num_connected,this->pmjs.size());
    for (uint32_t i = 0; i < this->pmjs.size(); i++)
        fprintf(file,"%g\n",this->pmjs[i]->error);
    fclose(file);
}

void PMJ_Data::write_errors_in_vtk (std::string filename)
{
    uint32_t np = this->pmjs.size();
    FILE *file = fopen(filename.c_str(),"w+");

    fprintf(file,"# vtk DataFile Version 4.1\n");
    fprintf(file,"vtk output\n");
    fprintf(file,"ASCII\n");
    fprintf(file,"DATASET POLYDATA\n");
    fprintf(file,"POINTS %u float\n",np);
    for (uint32_t i = 0; i < np; i++)
        fprintf(file,"%g %g %g\n",this->pmjs[i]->pos[0],this->pmjs[i]->pos[1],this->pmjs[i]->pos[2]);
        
    fprintf(file,"VERTICES %u %u\n",np,np*2);
    for (uint32_t i = 0; i < np; i++)
        fprintf(file,"1 %u\n",i);
    
    fprintf(file,"POINT_DATA %u\n",np);
    fprintf(file,"FIELD FieldData 1\n");
    fprintf(file,"LAT_error 1 %u float\n",np);
    for (uint32_t i = 0; i < np; i++)
        fprintf(file,"%g\n",this->pmjs[i]->error);
    
    fclose(file);
}

Surface_Data::Surface_Data (std::string filename)
{
    read_surface_from_file(filename);
}

Surface_Data::~Surface_Data ()
{
    
}

void Surface_Data::read_surface_from_file (std::string filename)
{
    vtkObject::GlobalWarningDisplayOff();

    vtkSmartPointer<vtkGenericDataObjectReader> reader = vtkSmartPointer<vtkGenericDataObjectReader>::New();
    reader->SetFileName(filename.c_str());
    reader->Update();

    if (reader->IsFilePolyData())
    {
        this->endo_data = reader->GetPolyDataOutput();
        this->dijkstra = vtkSmartPointer<vtkDijkstraGraphGeodesicPath>::New();
        this->dijkstra->SetInputData(this->endo_data);
        this->point_locator = vtkSmartPointer<vtkPointLocator>::New();
        this->point_locator->SetDataSet(this->endo_data);
        this->point_locator->BuildLocator();
    }
    else
    {
        fprintf(stderr,"[surface] Surface data must be Polydata!\n");
        exit(EXIT_FAILURE);
    }
}

std::vector<CloudPoint> Surface_Data::compute_geodesic_pathway (const double a[], const double b[], const bool use_first_node, const double m[])
{
    vtkIdType a_id = this->point_locator->FindClosestPoint(a);
    vtkIdType b_id = this->point_locator->FindClosestPoint(b);

    this->dijkstra->SetStartVertex(b_id);
    this->dijkstra->SetEndVertex(a_id);
    this->dijkstra->Update();

    std::vector<CloudPoint> geodesic_points;
    vtkSmartPointer<vtkIdList> ids = this->dijkstra->GetIdList();
    uint32_t num_points = ids->GetNumberOfIds();
    if (use_first_node)
    {
        for (uint32_t i = 0; i < num_points; i++)
        {
            double pos[3];
            uint32_t id = ids->GetId(i);
            this->endo_data->GetPoint(id,pos);

            CloudPoint cp(i,pos,0,false);
            geodesic_points.push_back(cp);
        }
    }
    else
    {
        CloudPoint cp(0,m,0,false);
        geodesic_points.push_back(cp);
        for (uint32_t i = 1; i < num_points; i++)
        {
            double pos[3];
            uint32_t id = ids->GetId(i);
            this->endo_data->GetPoint(id,pos);

            CloudPoint cp(i,pos,0,false);
            geodesic_points.push_back(cp);
        }
    }
    
    return geodesic_points;
}

void Surface_Data::print ()
{
    
}

void write_cloud_array_to_file (std::string filename, std::vector<CloudPoint*> arr)
{
    uint32_t np = arr.size();
    FILE *file = fopen(filename.c_str(),"w+");
    fprintf(file,"# vtk DataFile Version 4.1\n");
    fprintf(file,"vtk output\n");
    fprintf(file,"ASCII\n");
    fprintf(file,"DATASET POLYDATA\n");
    fprintf(file,"POINTS %u float\n",np);
    for (uint32_t i = 0; i < np; i++)
        fprintf(file,"%g %g %g\n",arr[i]->pos[0],arr[i]->pos[1],arr[i]->pos[2]);    
    fprintf(file,"VERTICES %u %u\n",np,np*2);
    for (uint32_t i = 0; i < np; i++)
        fprintf(file,"1 %u\n",i);
    fclose(file);
}