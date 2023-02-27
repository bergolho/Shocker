# ===========================================================================================================================
# This script generate the configuration files for all the Purkinje networks using different seed values.
# ---------------------------------------------------------------------------------------------------------------------------
# Author: Lucas Berg
# Last change: 02/08/2022
# ===========================================================================================================================

import benchmark as bc
import finsberg as fb
import canine as cn
import patient as pt
import oxford as ox
import sebastian as sb

def read_seeds_file (filename):
    seeds = []
    file = open(filename)
    for line in file:
        seeds.append(int(line))
    file.close()
    return seeds

def main ():
    seeds = read_seeds_file("seeds_100.txt")
    #seeds_LV = read_seeds_file("seeds_sebastian_LV_100.txt")
    #seeds_RV = read_seeds_file("seeds_sebastian_RV_100.txt")
    fb.write_finsberg_config_files(seeds)
    #cn.write_canine_config_files(seeds)
    #sb.write_sebastian_config_files(seeds_LV,seeds_RV)

if __name__ == "__main__":
    main()
