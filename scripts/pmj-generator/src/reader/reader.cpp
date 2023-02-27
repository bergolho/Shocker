#include "reader.h"

VTU_Reader::VTU_Reader (std::string filename)
{
    vtkSmartPointer<vtkXMLUnstructuredGridReader> reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
    reader->SetFileName(filename.c_str());
    reader->Update();

    vtkUnstructuredGrid *unstructured_grid = reader->GetOutput();

    uint32_t num_cells = unstructured_grid->GetNumberOfCells();
    uint32_t num_points = unstructured_grid->GetNumberOfPoints();

    for (uint32_t i = 0; i < num_points; i++)
    {
        double pos[3];
        unstructured_grid->GetPoint(i,pos);
        Point p(i,pos[0],pos[1],pos[2]);
        this->the_points.push_back(p);
    }

    for (uint32_t i = 0; i < num_cells; i++)
    {
        vtkCell *cell = unstructured_grid->GetCell(i);
        vtkHexahedron *hexahedron = dynamic_cast<vtkHexahedron*>(cell);

        uint32_t ids[8];
        double center[3];
        for (uint32_t j = 0; j < 8; j++)
        {
            ids[j] = hexahedron->GetPointIds()->GetId(j);
        }
        double dx = fabs(this->the_points[ ids[1] ].x - this->the_points[ ids[0] ].x);
        double dy = fabs(this->the_points[ ids[2] ].y - this->the_points[ ids[0] ].y);
        double dz = fabs(this->the_points[ ids[4] ].z - this->the_points[ ids[0] ].z);
        center[0] = this->the_points[ ids[0] ].x + (dx/2);
        center[1] = this->the_points[ ids[0] ].y + (dy/2);
        center[2] = this->the_points[ ids[0] ].z + (dz/2);
        
        Cell c(ids,center);
        this->the_cells.push_back(c);
    }

    std::string array_name = "Scalars_";
    vtkSmartPointer<vtkFloatArray> array = vtkFloatArray::SafeDownCast(unstructured_grid->GetCellData()->GetArray(array_name.c_str()));
    if(array)
    {
        for(uint32_t i = 0; i < num_cells; i++)
        {
            double value = array->GetValue(i);
            this->the_scalars.push_back(value);
        }
    }
}

void VTU_Reader::get_cells_center_positions(std::vector<Point> &out)
{
    for (uint32_t i = 0; i < this->the_cells.size(); i++)
    {
        double *center = this->the_cells[i].center;
        double lat = this->the_scalars[i];
        Point p(i,center[0],center[1],center[2],lat);
        out.push_back(p);
    }
}

void read_points_from_vtk (std::string filename, std::vector<Point> &arr)
{
    vtkSmartPointer<vtkGenericDataObjectReader> reader = vtkSmartPointer<vtkGenericDataObjectReader>::New();
    reader->SetFileName(filename.c_str());
    reader->Update();

    if(reader->IsFilePolyData())
    {
        vtkSmartPointer<vtkPolyData> input_polydata = reader->GetPolyDataOutput();
        uint32_t num_points = input_polydata->GetNumberOfPoints();

        for (uint32_t i = 0; i < num_points; i++)
        {
            double pos[3];
            input_polydata->GetPoint(i,pos);

            Point p(i,pos[0],pos[1],pos[2]);
            arr.push_back(p);
        }
        std::string array_name = "LAT";
        vtkSmartPointer<vtkFloatArray> array = vtkFloatArray::SafeDownCast(input_polydata->GetPointData()->GetArray(array_name.c_str()));
        if(array)
        {
            for(uint32_t i = 0; i < num_points; i++)
            {
                double value = array->GetValue(i);
                arr[i].value = value;
            }
        }
    }
    else
    {
        vtkSmartPointer<vtkUnstructuredGrid> input_grid = reader->GetUnstructuredGridOutput();
        uint32_t num_points = input_grid->GetNumberOfPoints();

        for (uint32_t i = 0; i < num_points; i++)
        {
            double pos[3];
            input_grid->GetPoint(i,pos);

            Point p(i,pos[0],pos[1],pos[2]);
            arr.push_back(p);
        }
        std::string array_name = "LAT";
        vtkSmartPointer<vtkFloatArray> array = vtkFloatArray::SafeDownCast(input_grid->GetPointData()->GetArray(array_name.c_str()));
        if(array)
        {
            for(uint32_t i = 0; i < num_points; i++)
            {
                double value = array->GetValue(i);
                arr[i].value = value;
            }
        }
    }
    
}

void read_pmjs (std::string filename, std::vector<PMJ> &pmjs)
{
    FILE *file = fopen(filename.c_str(),"r");
    uint32_t id, term_id;
    double pos[3], pk_lat, tiss_lat, delay, pk_dist, pk_cv;
    while (fscanf(file,"%u %lf %lf %lf %lf %lf %lf %u %lf %lf %lf",&term_id,&pos[0],&pos[1],&pos[2],&pk_lat,&tiss_lat,&delay,&id,&pk_dist,&pk_lat,&pk_cv) != EOF)
    {
        PMJ pmj(id,term_id,pos,pk_lat,tiss_lat,delay,pk_dist,pk_cv);
        pmjs.push_back(pmj);
    }
    fclose(file);
}

void get_closest_points (std::vector<Point> in, std::vector<Point> in2, std::vector<Point> &out, const uint32_t target_num_points)
{
    const uint32_t nc = in.size();
    std::vector<bool> taken;
    taken.assign(nc,false);

    for (uint32_t i = 0; i < in2.size(); i++)
    {
        //print_progress_bar(i,in2.size());
        
        uint32_t closest_id = 0;
        double closest_dist = __DBL_MAX__;
        Point cur_point = in2[i];
        for (uint32_t j = 0; j < nc; j += 100)
        {
            double dist = calc_norm(cur_point.x,cur_point.y,cur_point.z,in[j].x,in[j].y,in[j].z);
            if (dist < closest_dist)
            {
                closest_dist = dist;
                closest_id = j;
            }
        }
        taken[closest_id] = true;
    }

    std::vector<Point> points;
    for (uint32_t i = 0; i < nc; i++)
    {
        if (taken[i])
        {
            Point p(i,in[i].x,in[i].y,in[i].z,in[i].value);
            points.push_back(p);
        }
    }

    uint32_t np = points.size();
    uint32_t closest_np = target_num_points;
    
    taken.clear();
    taken.assign(np,false);

    uint32_t counter = 0;
    while (counter < closest_np)
    {
        uint32_t id = rand() % np;
        if (!taken[id])
        {
            out.push_back(points[id]);
            taken[id] = true;
            counter++;
        }
    }
}

void VTU_Reader::print ()
{
    printf("======= POINTS =======\n");
    for (uint32_t i = 0; i < this->the_points.size(); i++)
    {
        uint32_t id = this->the_points[i].id;
        double x = this->the_points[i].x;
        double y = this->the_points[i].y;
        double z = this->the_points[i].z;

        printf("[%u] = (%g %g %g)\n",id,x,y,z);
    }
 
    printf("======= CELLS =======\n");
    for (uint32_t i = 0; i < this->the_cells.size(); i++)
    {
        printf("[%u] -> center=(%g %g %g) {Scalar=%g}\n",i,this->the_cells[i].center[0],this->the_cells[i].center[1],this->the_cells[i].center[2],this->the_scalars[i]);
    }
}

std::vector<Point> filter_points (std::vector<Point> endo, const uint32_t target_num_points)
{
    uint32_t n = endo.size();
    uint32_t remaining_points = target_num_points;
    std::vector<Point> out;
    while (remaining_points > 0)
    {
        uint32_t selected_index = rand() % n;
        out.push_back(endo[selected_index]);
        remaining_points--;
    }
    return out;
}

std::vector<Point> filter_points_with_file (std::vector<Point> endo, std::vector<Point> pmjs, const uint32_t target_num_points)
{
    uint32_t n = endo.size();
    uint32_t remaining_points = target_num_points;
    std::vector<bool> taken;
    taken.assign(n,false);
    
    for (uint32_t i = 0; i < pmjs.size(); i++)
    {
        for (uint32_t j = 0; j < endo.size(); j++)
        {
            double dist = calc_norm(pmjs[i].x,pmjs[i].y,pmjs[i].z,\
                                    endo[j].x,endo[j].y,endo[j].z);
            if (dist < REGION_RADIUS)
            {
                taken[j] = true;
            }
        }
    }
    std::vector<Point> region_points;
    while (remaining_points > 0)
    {
        uint32_t id = rand() % n;
        if (taken[id])
        {
            region_points.push_back(endo[id]);
            remaining_points--;
        }
    }

    return region_points;
}

void concatenate_points (std::vector<Point> &out, std::vector<Point> in)
{
    for (uint32_t i = 0; i < in.size(); i++)
    {
        out.push_back(in[i]);
    }
}