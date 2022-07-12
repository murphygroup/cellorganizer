#!/bin/bash
# 
# 2022-06-23
# Copyright 2021-2022 Murphy Lab, CMU

# Determine the range of the number of vesicles from some collected data
# Define cellorganizer in your script before calling this script
if [ ! -v "cellorganizer" ]; then
    echo "Please set shell variable cellorganizer to the path of your CellOrganizer installation" >&2
    exit 1
fi
source "${cellorganizer}/module_if_available.sh" ; module_if_available load python36
read -r -d '' n_vesicles_python <<-'END_HEREDOC'
from math import *
'''Computed from data from these 102 samples ('generate_simulation_instances_min.20201006').'''
cytoplasm_volume_mean = 975.6973273003645
vesicle_volume_mean = 0.035302936430848024
'''Attempt to fill between 0 and cytoplasm_max_mean_fill_fraction of cytoplasm by setting `n_vesicles_min` and `n_vesicles_max` based on measurements from CellOrganizer-generated geometries. Intersecting objects will be discarded.'''
cytoplasm_max_mean_fill_fraction = 0.5
n_vesicles_max = ceil(cytoplasm_max_mean_fill_fraction * cytoplasm_volume_mean / vesicle_volume_mean)
n_vesicles_min = 1
#n_vesicles_min = ceil(n_vesicles_max / 10)
print(f'{n_vesicles_min} {n_vesicles_max}')
END_HEREDOC
n_vesicles=$(python3 - <<< "${n_vesicles_python}")
n_vesicles=(${n_vesicles})

