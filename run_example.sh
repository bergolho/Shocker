#!/bin/bash
# ===============================================================================================================================
# Author: Lucas Berg
# Last change: 20/07/2023
# This example demonstrate how to run the Shocker program and build a LVRV Purkinje network 
# for the Simplified mesh
# ===============================================================================================================================

SEED_NUMBER="1562072768"

if [ ! -f "./bin/Shocker" ]; then
    echo "[!] Compilation starts ..."
    ./recompile_project.sh
fi

# Generate the LV Purkinje network
echo "[!] Generating LV Purkinje network ..."
./bin/Shocker inputs/simplified_LV.ini
echo "[!] Copying the LV minimum network to the MergeNetwork folder"
cp ./outputs/simplified_LV_seed:${SEED_NUMBER}/minimum_network.vtk \
./scripts/merge-network/inputs/simplified_minimum_network_LV_seed:${SEED_NUMBER}.vtk

# Generate the RV Purkinje network
echo "[!] Generating RV Purkinje network ..."
./bin/Shocker inputs/simplified_RV.ini
echo "[!] Copying the RV minimum network to the MergeNetwork folder"
cp ./outputs/simplified_RV_seed:${SEED_NUMBER}/minimum_network.vtk \
./scripts/merge-network/inputs/simplified_minimum_network_RV_seed:${SEED_NUMBER}.vtk

# Merge the two networks in LVRV
if [ ! -f "./scripts/merge-network/bin/MergeNetwork" ]; then
    echo "[!] Compilation starts ..."
    cd ./scripts/merge-network
    ./recompile_project.sh
fi

# Generate the LVRV Purkinje network by merging the LV and RV networks using the His-bundle
echo "[!] Merging LVRV Purkinje network ..."
cd ./scripts/merge-network
./bin/MergeNetwork inputs/reference/simplified_his_um.vtk \
inputs/simplified_minimum_network_LV_seed:${SEED_NUMBER}.vtk \
inputs/simplified_minimum_network_RV_seed:${SEED_NUMBER}.vtk \
inputs/reference/simplified_root_positions_um.txt \
outputs/simplified_minimum_network_LVRV_seed:${SEED_NUMBER}.vtk

echo "[!] Merged LVRV Purkinje network saved at:> './scripts/merge-network/outputs'"