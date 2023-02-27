import numpy as np

def write_sebastian_LV_config_file (seed):
    filename = "config_files/LV/sebastian_LV_seed:%d.ini" % (seed)
    file = open(filename,"w")
    file.write('''[main]
#root_pos = [56136.2, 51915.3, 608]  ; Septum
root_pos = [39492.2, 50117.3, 93448]  ; Apex
seed = %d
start_diameter = 62.1468
l_d = 30000
his_offset = 0.0
use_initial_network = true
initial_network_filename = roots/sebastian_gold_LBB_um.txt
N_p = 20
N_a = 200

[save_network]
output_dir = outputs/sebastian_LV_seed:%d

[cloud_points]
phi = 0.01
rand_offset = 4
use_cloud_points = true
cloud_points_filename = clouds/inactives/sebastian_mesh_cloud_LV_um.vtk
surface_filename = clouds/endocardium/sebastian_endo_with_normals_LV_decimated_50perc_um.vtk 

[pmj]
use_pmj_location = true
pmj_location_filename = clouds/actives/sebastian_reference_pmjs_CV:1,9_LV_um.vtk
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

def write_sebastian_RV_config_file (seed):
    filename = "config_files/RV/sebastian_RV_seed:%d.ini" % (seed)
    file = open(filename,"w")
    file.write('''[main]
#root_pos = [29675.2, 73252.3, 12138]         ; Septum
root_pos = [41316.2, 99290.3, 51158]         ; Apex
seed = %d
start_diameter = 62.1468
l_d = 30000
his_offset = 0.0
use_initial_network = true
initial_network_filename = roots/sebastian_gold_RBB_um.txt
N_p = 20
N_a = 200

[save_network]
output_dir = outputs/sebastian_RV_seed:%d

[cloud_points]
phi = 0.01
rand_offset = 4
use_cloud_points = true
cloud_points_filename = clouds/inactives/sebastian_mesh_cloud_RV_um.vtk
surface_filename = clouds/endocardium/sebastian_endo_with_normals_RV_decimated_50perc_um.vtk

[pmj]
use_pmj_location = true
pmj_location_filename = clouds/actives/sebastian_reference_pmjs_CV:1,9_RV_um.vtk
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

def write_sebastian_config_files (seeds_LV, seeds_RV):
    for seed in seeds_LV:
        write_sebastian_LV_config_file(seed)
    for seed in seeds_RV:
        write_sebastian_RV_config_file(seed)
