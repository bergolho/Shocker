import numpy as np

def write_canine_LV_config_file (seed):
    filename = "config_files/LV/canine_LV_seed:%d.ini" % (seed)
    file = open(filename,"w")
    file.write('''[main]
root_pos = [39000.05, 60000.05, 13500]
seed = %d
start_diameter = 68.8608
l_d = 10000
his_offset = 0.0
N_p = 20
N_a = 120

[save_network]
output_dir = outputs/canine_LV_seed:%d

[cloud_points]
phi = 0.01
rand_offset = 4
use_cloud_points = true
cloud_points_filename = clouds/inactives/canine_mesh_cloud_LV_um.vtk
surface_filename = clouds/endocardium/canine_endo_smooth_with_normals_LV_decimated_90perc_um.vtk

[pmj]
use_pmj_location = true
pmj_location_filename = clouds/actives/canine_reference_pmjs_projected_CV:2_LV_um.vtk
pmj_connection_rate = 25
lat_error_tolerance = 2.0

[cost_function]
function_name = custom_function
beta = 1.0
alpha = 0.0
min_degrees_limit = 10
max_degrees_limit = 90
min_segment_length = 100
max_segment_length = 10000

    ''' % (seed,seed))
    file.close()

def write_canine_RV_config_file (seed):
    filename = "config_files/RV/canine_RV_seed:%d.ini" % (seed)
    file = open(filename,"w")
    file.write('''[main]
root_pos = [51500.05, 61750.05, 19750.05]
seed = %d
start_diameter = 68.8608
l_d = 10000
his_offset = 0.0
N_p = 20
N_a = 120

[save_network]
output_dir = outputs/canine_RV_seed:%d

[cloud_points]
phi = 0.01
rand_offset = 4
use_cloud_points = true
cloud_points_filename = clouds/inactives/canine_mesh_cloud_RV_um.vtk
surface_filename = clouds/endocardium/canine_endo_smooth_with_normals_RV_decimated_90perc_um.vtk

[pmj]
use_pmj_location = true
pmj_location_filename = clouds/actives/canine_reference_pmjs_projected_CV:2_RV_um.vtk
pmj_connection_rate = 25
lat_error_tolerance = 2.0

[cost_function]
function_name = custom_function
beta = 1.0
alpha = 0.0
min_degrees_limit = 10
max_degrees_limit = 90
min_segment_length = 100
max_segment_length = 10000

    ''' % (seed,seed))
    file.close()   

def write_canine_config_files (seeds):
    for seed in seeds:
        write_canine_LV_config_file(seed)
        write_canine_RV_config_file(seed)
