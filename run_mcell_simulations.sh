#!/bin/bash
# Run MCell simulations locally or on the cluster

usage()
{
    cat << EOF
usage: run_mcell_simulations.sh --help | [--dry_run] [--mode <mode>] [--cluster_partition <partition>] [--cluster_cpus <ncpus>] [--cluster_memory <MB>] [--cluster_exclusive] [--cluster_jobs <njobs>] [--seed_offset <off>] [--verbose] [--keep_cluster_scripts] [--ignore_errors] [--keep_existing_results] directory_of_simulations ...
EOF
}

#directory_of_simulations=()
declare -a directory_of_simulations
dry_run=0
mode="local"
n_seeds=1
cluster_partition=""
cluster_cpus=""
cluster_memory=""
cluster_jobs=1
cluster_exclusive=0
seed_offset=0
verbose=0
keep_cluster_scripts=0
ignore_errors=0
keep_existing_results=0

option_pattern='^--[[A-Za-z0-9_\-]\+$'
options_pattern='^--\(help\|dry_run\|mode\|n_seeds\|cluster_partition\|cluster_cpus\|cluster_memory\|cluster_exclusive\|cluster_jobs\|seed_offset\|verbose\|keep_cluster_scripts\|ignore_errors\|keep_existing_results\)$'
mode_pattern='^\(local\|slurm\)$'

integer_pattern='^[0-9]\+$'
cluster_partition_pattern='^[A-Za-z0-9][A-Za-z0-9_]*$'

timestamp_cmd="date +%Y%m%d%H%M%S"
timestamp_cmd_long="date +\"%Y-%m-%d %H:%M:%S\""
timestamp_cmd_long_quoted="date +\\\"%Y-%m-%d %H:%M:%S\\\""

concatenate_array()
{
    args2=("$@")
    numargs2="${#args2[@]}"
    separator="${args2[${numargs2}-1]}"
    # exit 1
    for i in $(seq 0 ${numargs2-2}); do
        if [ $i -eq "0" ]; then
            echo -n "${args2[${i}]}"
        else
            echo -n "$separator${args2[${i}]}"
        fi
    done
    echo " "
}

execute_command()
{
    if [ "$dry_run" = 1 ]; then
        echo "Dry run; would have executed '$@'"
    else
        echo "Executing '$@'"
        eval $@
    fi
}


if [ "$#" -eq "0" ]; then
    usage
    exit 1
fi

while (( "$#" )); do
    if [ `expr match "$1" "$option_pattern"` -gt "0" ]; then
        option=`echo "$1" | cut -c3-`
        shift
        
        case "$option" in
            help)
                usage
                ;;
                
            dry_run)
                dry_run=1
                ;;
                
            mode)
                if [ "$#" -lt 1 ]; then
                    echo "Option --${option} requires one argument" >&2
                    exit 1
                fi
                if [ `expr length "\`expr match \"$1\" \"$mode_pattern\"\`"` -eq "0" ]; then
                    echo "Invalid mode ${1}" >&2
                    exit 1
                fi
                mode="$1"
                shift
                ;;
                
            n_seeds)
                if [ "$#" -lt 1 ]; then
                    echo "Option --${option} requires one argument" >&2
                    exit 1
                fi
                if [ `expr length "\`expr match \"$1\" \"$integer_pattern\"\`"` -eq "0" ]; then
                    echo "Invalid cluster_cpus ${1}" >&2
                    exit 1
                fi
                n_seeds="$1"
                shift
                ;;
            
            cluster_partition)
                if [ "$#" -lt 1 ]; then
                    echo "Option --${option} requires one argument" >&2
                    exit 1
                fi
                if [ `expr length "\`expr match \"$1\" \"$cluster_partition_pattern\"\`"` -eq "0" ]; then
                    echo "Invalid cluster_partition ${1}" >&2
                    exit 1
                fi
                cluster_partition="$1"
                shift
                ;;
            
            cluster_cpus)
                if [ "$#" -lt 1 ]; then
                    echo "Option --${option} requires one argument" >&2
                    exit 1
                fi
                if [ `expr length "\`expr match \"$1\" \"$integer_pattern\"\`"` -eq "0" ]; then
                    echo "Invalid cluster_cpus ${1}" >&2
                    exit 1
                fi
                cluster_cpus="$1"
                shift
                ;;
            
            cluster_memory)
                if [ "$#" -lt 1 ]; then
                    echo "Option --${option} requires one argument" >&2
                    exit 1
                fi
                if [ `expr length "\`expr match \"$1\" \"$integer_pattern\"\`"` -eq "0" ]; then
                    echo "Invalid cluster_memory ${1}" >&2
                    exit 1
                fi
                cluster_memory="$1"
                shift
                ;;
            
            cluster_exclusive)
                cluster_exclusive=1
                ;;
            
            matlab_name)
                if [ "$#" -lt 1 ]; then
                    echo "Option --${option} requires one argument" >&2
                    exit 1
                fi
                if [ `expr length "\`expr match \"$1\" \"$integer_pattern\"\`"` -eq "0" ]; then
                    echo "Invalid matlab_name ${1}" >&2
                    exit 1
                fi
                matlab_name="$1"
                shift
                ;;
            
            cluster_jobs)
                if [ "$#" -lt 1 ]; then
                    echo "Option --${option} requires one argument" >&2
                    exit 1
                fi
                if [ `expr length "\`expr match \"$1\" \"$integer_pattern\"\`"` -eq "0" ]; then
                    echo "Invalid cluster_jobs ${1}" >&2
                    exit 1
                fi
                cluster_jobs="$1"
                shift
                ;;
            
            seed_offset)
                if [ "$#" -lt 1 ]; then
                    echo "Option --${option} requires one argument" >&2
                    exit 1
                fi
                if [ `expr length "\`expr match \"$1\" \"$integer_pattern\"\`"` -eq "0" ]; then
                    echo "Invalid seed_offset ${1}" >&2
                    exit 1
                fi
                seed_offset="$1"
                shift
                ;;
            
            verbose)
                verbose=1
                ;;
            
            keep_cluster_scripts)
                keep_cluster_scripts=1
                ;;
            
            ignore_errors)
                ignore_errors=1
                ;;
            
            keep_existing_results)
                keep_existing_results=1
                ;;
            
            *)
                echo "Invalid option --${option}" >&2
                exit 1
                ;;
        esac
        
    else
        shopt -s nullglob
        directory_of_simulations+=($1)
        shopt -u nullglob
        shift
        
    #else
        #echo "Too many arguments" >&2
        #exit 1
    fi
done

mcell_main_scripts=()
for i in $(seq 0 $(expr ${#directory_of_simulations[@]} - 1)); do
    mcell_main_scripts+=($(find "${directory_of_simulations[i]}" -regex '.*/cell\(\.main\)?\.\(mdl\|mcell\)' | sort))
done
n_mcell_main_scripts="${#mcell_main_scripts[@]}"

if (( $keep_existing_results )) ; then
    # Do not run scripts with react_data or viz_data subdirectories
    shopt -s globstar
    for i in $(seq 0 $(expr $n_mcell_main_scripts - 1)); do
        mcell_script="${mcell_main_scripts[${i}]}"
        mcell_script_dir="'${mcell_script}'/cell\(\.main\)?\.\(mdl\|mcell\)//"
        mcell_script_react_data_pattern="${mcell_script_dir}react_data/**/*.dat"
        mcell_script_viz_data_pattern="${mcell_script_dir}viz_data/**/*.dat"
        
        shopt -s nullglob
        mcell_script_react_data_files=($mcell_script_react_data_pattern)
        mcell_script_viz_data_files=($mcell_script_viz_data_pattern)
        shopt -u nullglob
        
        # FIXME: This might break with spaces in filenames
        for f in "${mcell_script_react_data_files[@]}" ; do
            unset mcell_main_scripts[$i]
            break
        done
        for f in "${mcell_script_viz_data_files[@]}" ; do
            unset mcell_main_scripts[$i]
            break
        done
    done
    mcell_main_scripts=("${mcell_main_scripts[@]}")
fi
n_mcell_main_scripts="${#mcell_main_scripts[@]}"
case "$mode" in
    slurm)
    slurm_output_base='${SLURM_SUBMIT_DIR}/slurm_job.${SLURM_JOBID}.${HOST}.${USER}'
    slurm_stdout="${slurm_output_base}.stdout"
    slurm_stderr="${slurm_output_base}.stderr"
    slurm_command_base="/usr/bin/sbatch"
    if [ `expr length "$cluster_partition"` -gt "0" ]; then
        slurm_command_base="${slurm_command_base} --partition=${cluster_partition}"
    fi
    if [ `expr length "$cluster_cpus"` -gt "0" ]; then
        slurm_command_base="${slurm_command_base} --cpus-per-task=${cluster_cpus}"
    fi
    if [ `expr length "$cluster_memory"` -gt "0" ]; then
        slurm_command_base="${slurm_command_base} --mem=${cluster_memory}"
    fi
    if [ "$cluster_exclusive" = 1 ]; then
        slurm_command_base="${slurm_command_base} --exclusive"
    fi
esac

original_dir=$(pwd)
tempdir="${HOME}/tmp"
if [ ! -e "${tempdir}" ]; then
    mkdir "${tempdir}"
fi
let n_mcell_commands="$n_mcell_main_scripts * $n_seeds"
case "$mode" in
    local)
        let n_mcell_commands_per_job=n_mcell_commands
        ;;
    
    slurm)
        let n_mcell_commands_per_job="( n_mcell_commands + cluster_jobs - 1 ) / cluster_jobs"
        ;;
esac
mcell_commands=()
n_mcell_commands_done=0
job_count=0
tempfilenames=()
tempfilename=""
execute_command_suffix=""
if [ "$ignore_errors" = 1 ]; then
    execute_command_suffix=" || :"
fi
for i in $(seq 0 $(expr $n_mcell_main_scripts - 1)); do
    mcell_script="${mcell_main_scripts[${i}]}"
    for j in $(seq 0 $(expr $n_seeds - 1)); do
        if (( $n_mcell_commands_done % $n_mcell_commands_per_job == 0 )) ; then
            tempfilename=$(mktemp "--tmpdir=${tempdir}")
            tempfilenames+=($tempfilename)
            : > "${tempfilename}"
            cat <<- EOF >> "${tempfilename}"
#!/bin/bash
HOST=\`hostname -s\`
USER=\`whoami\`
EOF
            let job_count++
        fi

        mcell_command_index=$n_mcell_commands_done
        cmd_ind=$mcell_command_index
        cat <<- EOF >> "${tempfilename}"

MCELL_SCRIPT="${mcell_script}"
echo "Running MCell on \${MCELL_SCRIPT}"
MCELL_SCRIPT_DIR="\$(dirname "\${MCELL_SCRIPT}")"
MCELL_SEED=$(expr $cmd_ind + 1 + ${seed_offset})
MCELL_SCRIPT_STDOUT="\${MCELL_SCRIPT}.\${MCELL_SEED}.stdout"
MCELL_SCRIPT_STDERR="\${MCELL_SCRIPT}.\${MCELL_SEED}.stderr"
echo > "\$MCELL_SCRIPT_STDOUT"
function echo2 {
    echo \$@ >> "\$MCELL_SCRIPT_STDOUT"
}
MCELL_DRY_RUN=0
# MCELL_DRY_RUN=1
execute_command()
{
    if (( \$MCELL_DRY_RUN )) ; then
        echo2 "Dry run; would have executed 'eval \$@'"
    else
        echo2 eval \$@
        eval \$@
    fi
}
execute_command_with_optional_ignore_error()
{
    execute_command \$@ $execute_command_suffix
}
execute_command cd "\$MCELL_SCRIPT_DIR"
echo2
echo2 pwd=\`pwd\`
echo2 SLURM_SUBMIT_DIR=\$SLURM_SUBMIT_DIR
echo2 SLURM_JOBID=\$SLURM_JOBID
echo2 HOST=\$HOST
echo2 USER=\$USER
echo2 MCELL_SCRIPT=\$MCELL_SCRIPT
echo2 MCELL_SCRIPT_DIR=\$MCELL_SCRIPT_DIR
echo2 MCELL_SEED=\$MCELL_SEED
echo2 MCELL_SCRIPT_STDOUT=\$MCELL_SCRIPT_STDOUT
echo2 MCELL_SCRIPT_STDERR=\$MCELL_SCRIPT_STDERR
echo2
echo2 Started at \`$timestamp_cmd_long\`
echo2
echo2
execute_command_with_optional_ignore_error mcell -seed \$MCELL_SEED "\$MCELL_SCRIPT" 1>> "\$MCELL_SCRIPT_STDOUT" 2> "\$MCELL_SCRIPT_STDERR"
echo2
echo2
echo2 Finished at \`$timestamp_cmd_long\`
execute_command cd "${original_dir}"
EOF
        
        chmod a+rx "${tempfilename}"
        let n_mcell_commands_done++
    done
done

job_name_base="mcell_`$timestamp_cmd`"
for tempfilename in ${tempfilenames[@]} ; do
    echo
    if [ "$verbose" = 1 ]; then
        cat "${tempfilename}"
    fi
    echo
    
    case "$mode" in
        local)
            execute_command bash "${tempfilename}"
            ;;
        
        slurm)
            job_name="${job_name_base}_${job_count}"
            let job_count++
            slurm_command="${slurm_command_base} --job-name=\"${job_name}\" \"${tempfilename}\""
            execute_command "${slurm_command}"
            ;;
    esac
    
    if [ "$keep_cluster_scripts" = 0 ]; then
        rm "${tempfilename}"
    fi
done

