#!/bin/bash
# 
# 2021-08-04
# Copyright 2021-2022 Murphy Lab, CMU

###########################################################################
# SET THESE VARIABLES IN YOUR SCRIPT

## Arguments and brief documentation can be found by running `python cellorganizer/applications/generate_and_simulate.py --help`
# The directory into which you have installed CellOrganizer. Replace the relative path or replace everthing after `=` with an absolute path.
#cellorganizer=$(realpath '../cellorganizer')
# The name of the directory in which you have placed your `.mdl` files and the prefix for those files. Put this directory in the `data` subdirectory of your CellOrganizer installation.
#reaction_network='{data}/CBExMinScaled3EN20min'
# With no vesicles (`synthesis == 'framework'`)
#reaction_network='{data}/CBExMinScaled3_20min'
#framework_model='{models}/3D/spharm/lamp2.mat'
#vesicle_model='{models}/3D/tfr.mat'
#n_images_to_synthesize=100
#cluster_mode='slurm'

#synthesis="framework"
#downsampling="0.5"
#vcml_relative_downsampling="0.5"
#simulation_end_time="4000"

###########################################################################
# DO NOT MODIFY THIS BLOCK

declare -A required_vars
required_vars['cellorganizer']='the path of your CellOrganizer installation'
#required_vars['reaction_network_pattern']='either a path of a single VCML file or a glob pattern for a single collection of MCell MDL files'
required_vars['reaction_network_pattern']='a glob pattern for a single collection of MCell MDL files'
#required_vars['output_dir']='the path of the directory into which output should be written'
for required_var in "${!required_vars[@]}" ; do
    if [ ! -v "$required_var" ]; then
        echo "Please set shell variable $required_var to ${required_vars[$required_var]}" >&2
        exit 1
    fi
done

if [ ! -v "output_dir" ]; then
    if [ -n "$SLURM_JOBID" ]; then
        calling_script=$(scontrol show job "$SLURM_JOBID" | grep -E '^ *Command=' | sed -e 's/^ *Command=//')
    else
        calling_script=$(realpath "$0")
    fi
    #output_dir="${calling_script%.sh}"
    output_dir="${calling_script}.out"
fi

co_apps=$(realpath "${cellorganizer}/applications")
#co_data=$(realpath "${cellorganizer}/data")
#co_models=$(realpath "${cellorganizer}/models")
#https://stackoverflow.com/questions/56962129/how-to-get-original-location-of-script-used-for-slurm-job
cd "$co_apps"
# Uncomment either the MDL or the VCML line to use MCell or Virtual Cell
#reaction_network_pattern="${co_data}/${reaction_network}/${reaction_network}.*.mdl"
#reaction_network_pattern="${co_data}/${reaction_network}.vcml"
timestamp_pretty=$(date -u "+%Y-%m-%d %H:%M:%S %Z")
timestamp="${timestamp_pretty//[: -]/}"
echo "Beginning at $timestamp_pretty"


# debug settings
#n_images_to_synthesize=2 # debug
#n_images_to_synthesize=4 # debug
#n_images_to_synthesize=12 # debug
#cluster_mode='local' # debug
#simulation_end_time="40" # debug
#dry_run=1 # debug
#run_simulations=0 # debug
#run_analysis=0 # debug
#overwrite=1 # debug
#overwrite_analysis=1 # debug


#echo "\$@=$@"

# Generate geometries and run simulations

if [ ! -a "$output_dir" ]; then
    mkdir -p "$output_dir"
fi

function echo_if_def_true()
{
    local arg2
    if [[ "$#" > 1 ]]; then
        arg2="$2"
    else
        arg2="$1"
    fi
    if [[ -v "$1" && "${!1}" != "0" ]]; then
        #echo "--$arg2"
        echo -n "--$arg2"
    fi
}
#unset asdf ; echo "asdf='$asdf'" ; echo_if_def_true asdf ; asdf='0' ; echo "asdf='$asdf'" ; echo_if_def_true asdf ; asdf='1' ; echo "asdf='$asdf'" ; echo_if_def_true asdf; asdf='1' ; echo "asdf='$asdf'" ; echo_if_def_true asdf qwerty ; exit 1 # Debug
function append_option_if_def_true()
{
    # Argument 1 is array name to which option should be appended
    # Argument 2 is option name to be appended
    local arg1_name arg2_name arg2
    arg1_name="$1"
    arg2_name="$2"
    # Sanitize
    arg1_name="${arg1_name//[^A-Za-z0-9_]/}"
    arg2_name="${arg2_name//[^A-Za-z0-9_]/}"
    arg2="$arg2_name"
    if [[ -v "$arg2_name" && "${!arg2_name}" != "0" ]]; then
        eval "${arg1_name}+=(\"--\${arg2_name}\")"
    fi
}

function echo_value_if_def()
{
    # Argument 1 is option name to be appended
    # Argument 2 (optional) is variable name to check instead
    local arg2 arg2_deref
    if [[ "$#" > 1 ]]; then
        arg2="$2"
    else
        arg2="$1"
    fi
    if [[ -v "$1" && "${!1}" != "0" ]]; then
        arg2_deref="${!2}"
        echo -n "--$1 "
        printf "%q" "$arg2_deref"
    fi
}
#unset asdf ; echo "asdf='$asdf'" ; echo_value_if_def asdf ; asdf='123' ; echo "asdf='$asdf'" ; echo_value_if_def asdf ; asdf='123' ; echo "asdf='$asdf'" ; echo_value_if_def asdf qwerty ; exit 1  # Debug
#echo "synthesis='$synthesis'" ; echo_value_if_def synthesis ; echo_value_if_def synthesis synthesis_renamed ; exit 1  # Debug
function append_option_with_value_if_def()
{
    # Argument 1 is array name to which option should be appended
    # Argument 2 is option name to be appended
    local arg1_name arg2_name arg2
    arg1_name="$1"
    arg2_name="$2"
    # Sanitize
    arg1_name="${arg1_name//[^A-Za-z0-9_]/}"
    arg2_name="${arg2_name//[^A-Za-z0-9_]/}"
    arg2="$arg2_name[@]"
    if [[ -v "$arg2_name" ]]; then
        arg2=("${!arg2}")
        eval "${arg1_name}+=(\"--\${arg2_name}\")"
        eval "${arg1_name}+=(\"\${arg2[@]}\")"
    fi
}

function not_zero()
{
    [[ ! -v "$1" || "${!1}" -ne "0" ]]
}
#unset asdf ; echo "asdf='$asdf'" ; not_zero asdf ; echo "$?" ; asdf='0' ; echo "asdf='$asdf'" ; not_zero asdf ; echo "$?" ; asdf='123' ; echo "asdf='$asdf'" ; not_zero asdf ; echo "$?" ; exit 1  # Debug

function exit_code_to_boolean()
{
    [[ "$1" -ne "0" ]]
    echo "$?"
}

#echo '@@@@@@@@@@@@ DEBUG exiting early' ; exit 1 # Debug

generate_and_simulate_successful=0
if not_zero run_simulations; then
    args=()
    args+=("$reaction_network_pattern")
    args+=("$output_dir")
    append_option_if_def_true args dry_run
    append_option_if_def_true args overwrite overwrite_simulations
    append_option_with_value_if_def args cellorganizer
    append_option_with_value_if_def args matlab_setup
    append_option_with_value_if_def args cluster_partition
    append_option_with_value_if_def args generation_cluster_jobs
    append_option_with_value_if_def args generation_cluster_memory
    append_option_if_def_true args generation_cluster_exclusive
    append_option_with_value_if_def args simulation_cluster_jobs
    append_option_with_value_if_def args simulation_cluster_cpus
    append_option_with_value_if_def args simulation_cluster_memory
    append_option_with_value_if_def args framework_model
    append_option_with_value_if_def args vesicle_models
    append_option_with_value_if_def args synthesis
    append_option_if_def_true args output_image
    append_option_with_value_if_def args n_images_to_synthesize
    append_option_with_value_if_def args n_vesicles
    append_option_with_value_if_def args vesicle_volume_scale
    append_option_with_value_if_def args downsampling
    append_option_with_value_if_def args vcml_relative_downsampling
    append_option_with_value_if_def args framework_min_clearance
    append_option_with_value_if_def args framework_clearance_n_max_filter_rounds
    #append_option_if_def_true args remove_intersecting_mesh_objects
    append_option_with_value_if_def args intersecting_mesh_object_policy
    append_option_with_value_if_def args base_seed
    append_option_with_value_if_def args n_vesicles_seed_offset
    append_option_with_value_if_def args translations
    append_option_with_value_if_def args simulation_end_time
    append_option_with_value_if_def args simulation_default_time_step
    append_option_with_value_if_def args simulation_max_time_step
    append_option_with_value_if_def args simulation_output_time_step
    append_option_with_value_if_def args simulation_absolute_tolerance
    append_option_with_value_if_def args simulation_relative_tolerance
    append_option_with_value_if_def args simulation_interaction_radius
    append_option_with_value_if_def args simulation_seed_offset
    (
        set -o pipefail
        python3.6 -m IPython --pdb -- generate_and_simulate.py "${args[@]}" \
            2>&1 | tee generate_and_simulate.py.${timestamp}.log
    )
    generate_and_simulate_exit_code="$?"
    generate_and_simulate_successful=$(exit_code_to_boolean "$generate_and_simulate_exit_code")
    if (( ! generate_and_simulate_successful )); then
        echo 'generate_and_simulate.py failed' ; exit 1
    fi
fi

# Analyze results, create figures and tables
if (( generate_and_simulate_successful && run_analysis )); then
    cd "$co_apps"
    args=()
    args+=("$reaction_network_pattern")
    args+=("$output_dir")
    append_option_if_def_true args dry_run
    append_option_if_def_true args overwrite overwrite_analysis
    append_option_with_value_if_def args cellorganizer
    append_option_with_value_if_def args synthesis
    append_option_with_value_if_def args downsampling
    append_option_with_value_if_def args vcml_relative_downsampling
    append_option_with_value_if_def args simulation_end_time
    (
        set -o pipefail
        python3.6 -m IPython --pdb -- generate_and_simulate_analysis.py "${args[@]}" \
            2>&1 | tee generate_and_simulate_analysis.py.${timestamp}.log
    )
    generate_and_simulate_analysis_exit_code="$?"
    generate_and_simulate_analysis_successful=$(exit_code_to_boolean "$generate_and_simulate_analysis_exit_code")
    if (( ! generate_and_simulate_analysis_successful )); then
        echo 'generate_and_simulate_analysis.py failed' ; exit 1
    fi
fi
