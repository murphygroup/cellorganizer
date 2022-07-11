#!/bin/bash
# Run a Matlab script locally or on the cluster

args=("$@")

usage()
{
    cat << EOF
usage: run_matlab_script.sh --help | [--dry_run] [--mode <mode>] [--cluster_partition <partition>] [--cluster_cpus <ncpus>] [--cluster_memory <MB>] [--cluster_exclusive] [--cluster_jobs <njobs>] [--job_name <name>] [--matlab_location <location>] [--matlab_setup <script>] matlab_script
EOF
}

dry_run=0
mode="local"
matlab_script=""
job_name="matlab script"
cluster_partition=""
cluster_cpus=""
cluster_memory=""
cluster_jobs=1
cluster_exclusive=0
matlab_location="matlab"
matlab_setup=""
# matlab_setup="module load matlab-9.7"

option_pattern='^--[[A-Za-z0-9_\-]\+$'
options_pattern='^--\(help\|dry_run\|mode\|cluster_partition\|cluster_cpus\|cluster_memory\|cluster_exclusive\|cluster_jobs\|job_name\|\|matlab_location\|matlab_setup\)$'
mode_pattern='^\(local\|slurm\)$'
application_pattern='^[A-Za-z][A-Za-z0-9_\-\.]*$'

integer_pattern='^[0-9]\+$'
cluster_partition_pattern='^[A-Za-z0-9][A-Za-z0-9_]*$'

concatenate_array()
{
    args2=("$@")
    numargs2="${#args2[@]}"
    separator="${args2[${numargs2}-1]}"
    # exit 1
    for i in $(seq 0 ${numargs2-2}); do
        if [ $i -eq "0" ]; then
            /bin/echo -n "${args2[${i}]}"
        else
            /bin/echo -n "$separator${args2[${i}]}"
        fi
    done
    /bin/echo " "
}

execute_command()
{
    if [ "$dry_run" = 1 ]; then
        echo "Dry run; would have executed '$@'"
    else
        echo "Executing '$@'"
        #$@
        eval $@
    fi
}


if [ "$#" -eq "0" ]; then
    usage
    exit 1
fi

while (( "$#" )); do
    if [ ! -z "`expr match \"$1\" \"$options_pattern\"`" ]; then
        option=`/bin/echo "$1" | cut -c3-`
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
                    /bin/echo "Option --${option} requires one argument" >&2
                    exit 1
                fi
                if [ `expr length "\`expr match \"$1\" \"$mode_pattern\"\`"` -eq "0" ]; then
                    /bin/echo "Invalid mode ${arg0}" >&2
                    exit 1
                fi
                mode="$1"
                shift
                ;;
                
            cluster_partition)
                if [ "$#" -lt 1 ]; then
                    /bin/echo "Option --${option} requires one argument" >&2
                    exit 1
                fi
                if [ `expr length "\`expr match \"$1\" \"$cluster_partition_pattern\"\`"` -eq "0" ]; then
                    /bin/echo "Invalid cluster_partition ${arg0}" >&2
                    exit 1
                fi
                cluster_partition="$1"
                shift
                ;;
            
            cluster_cpus)
                if [ "$#" -lt 1 ]; then
                    /bin/echo "Option --${option} requires one argument" >&2
                    exit 1
                fi
                if [ `expr length "\`expr match \"$1\" \"$integer_pattern\"\`"` -eq "0" ]; then
                    /bin/echo "Invalid cluster_cpus ${arg0}" >&2
                    exit 1
                fi
                cluster_cpus="$1"
                shift
                ;;
            
            cluster_memory)
                if [ "$#" -lt 1 ]; then
                    /bin/echo "Option --${option} requires one argument" >&2
                    exit 1
                fi
                if [ `expr length "\`expr match \"$1\" \"$integer_pattern\"\`"` -eq "0" ]; then
                    /bin/echo "Invalid cluster_memory ${arg0}" >&2
                    exit 1
                fi
                cluster_memory="$1"
                shift
                ;;
            
            cluster_exclusive)
                cluster_exclusive=1
                ;;
            
            matlab_location)
                if [ "$#" -lt 1 ]; then
                    /bin/echo "Option --${option} requires one argument" >&2
                    exit 1
                fi
                matlab_location="$1"
                shift
                ;;
            
            job_name)
                if [ "$#" -lt 1 ]; then
                    /bin/echo "Option --${option} requires one argument" >&2
                    exit 1
                fi
                job_name="$1"
                shift
                ;;
            
            matlab_setup)
                if [ "$#" -lt 1 ]; then
                    /bin/echo "Option --${option} requires one argument" >&2
                    exit 1
                fi
                matlab_setup="$1"
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
            
            *)
                /bin/echo "Invalid option --${option}" >&2
                exit 1
                ;;
        esac
        
    elif [ `expr length "$matlab_script"` -eq "0" ]; then
        matlab_script="$1"
        shift
        
    else
        /bin/echo "Too many arguments" >&2
        exit 1
    fi
done


matlab_bash_command=""
if [ -n "$matlab_setup" ]; then
    matlab_bash_command+="${matlab_setup} ; "
fi
matlab_bash_command+="time ${matlab_location} -nosplash -nodesktop -singleCompThread -r \"try; setup; ${matlab_script}; catch err; disp(getReport(err, 'extended')); end; exit\""

case "$mode" in
    local)
        execute_command "${matlab_bash_command}"
        ;;
        
    slurm)
        slurm_output_base='${SLURM_SUBMIT_DIR}/slurm_job.${SLURM_JOBID}.${HOST}.${USER}'
        slurm_stdout="${slurm_output_base}.stdout"
        slurm_stderr="${slurm_output_base}.stderr"
        slurm_command="/usr/bin/sbatch"
        if [ `expr length "$cluster_partition"` -gt "0" ]; then
            slurm_command="${slurm_command} --partition=${cluster_partition}"
        fi
        if [ `expr length "$cluster_cpus"` -gt "0" ]; then
            slurm_command="${slurm_command} --cpus-per-task=${cluster_cpus}"
        fi
        if [ `expr length "$cluster_memory"` -gt "0" ]; then
            slurm_command="${slurm_command} --mem=${cluster_memory}"
        fi
        if [ "$cluster_exclusive" = 1 ]; then
            slurm_command="${slurm_command} --exclusive"
        fi
        
        tempdir="${HOME}/tmp"
        if [ ! -e "${tempdir}" ]; then
            mkdir "${tempdir}"
        fi
        tempfilename=$(mktemp "--tmpdir=${tempdir}")
        cat <<- EOF > "${tempfilename}"
#!/bin/sh
HOST=\`hostname -s\`
USER=\`whoami\`
pwd
echo \$SLURM_SUBMIT_DIR
echo \$SLURM_JOBID
echo \$HOST
echo \$USER
${matlab_bash_command} 1>"${slurm_stdout}" 2>"${slurm_stderr}"
EOF
        #cat "${tempfilename}"
        slurm_command="${slurm_command} --job-name=\"${job_name}\" \"${tempfilename}\""
        for i in $(seq 1 ${cluster_jobs}); do
            if (( i > 1 )); then
                sleep 1
            fi
            execute_command "${slurm_command}"
        done
        rm "${tempfilename}"
        ;;
        
    *)
        /bin/echo "Invalid mode $mode" >&2
        exit 1
        ;;
esac

