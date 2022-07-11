#!/bin/bash
# 
# 2021-08-04
# Copyright 2021-2022 Murphy Lab, CMU

###########################################################################
# FEEL FREE TO MODIFY THE VARIABLES IN THIS BLOCK

gas_args=()
gas_args+=
# The directory into which you have installed CellOrganizer. Replace the relative path or replace everthing after `=` with an absolute path.
cellorganizer=$(realpath '../cellorganizer')
# The name of the directory in which you have placed your `.mdl` files and the prefix for those files. Put this directory in the `data` subdirectory of your CellOrganizer installation.
#reaction_network_pattern='{data}/CBExMinScaled3EN20min.vcml'
# With no vesicles (`synthesis == 'framework'`)
reaction_network_pattern='{data}/CBExMinScaled3_05min/CBExMinScaled3_05min.*.mdl'
# New model is now the default
#framework_model='{models}/3D/spharm/lamp2.demo3D52.mat'
vesicle_model='{models}/3D/tfr.mat'
# Not using vesicles for now
unset vesicle_model # Debug
n_images_to_synthesize=5
#cluster_mode='slurm'
cluster_mode='local'
generation_cluster_jobs=2

synthesis='framework'
downsampling=0.5
vcml_relative_downsampling=0.5
framework_min_clearance=0.1

# Determine the range of the number of vesicles from some collected data
. applications/compute_n_vesicles_20220623.sh
echo "n_vesicles='${n_vesicles[@]}'" # Debug


simulation_end_time=4000

#run_simulations=0
run_simulations=1
#run_analysis=0
run_analysis=1


# Debug
n_images_to_synthesize=1
overwrite_simulations=1

###########################################################################
# DO NOT MODIFY THIS BLOCK

. "${cellorganizer}/applications/generate_and_simulate.sh"
