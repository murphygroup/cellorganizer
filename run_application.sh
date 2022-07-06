#!/bin/bash
# Run a application locally or on the cluster

args=("$@")

usage()
{
    cat << EOF
usage: run_application.sh --help | [--dry_run] [--mode <mode>] [--cluster_partition <partition>] [--cluster_cpus <ncpus>] [--cluster_memory <MB>] [--cluster_exclusive] [--cluster_jobs <njobs>] [--matlab_name <name>] [--matlab_setup <command>] [--do_not_delete_temp] application_to_run
EOF
}

dry_run=0
mode="local"
application_to_run=""
cluster_partition=""
cluster_cpus=""
cluster_memory=""
cluster_jobs=1
cluster_exclusive=0
matlab_name="matlab"
# matlab_name="matlab_8.6"
matlab_setup=":"
# matlab_setup="module load matlab-9.7"
delete_temp=1

option_pattern='^--[[A-Za-z0-9_\-]\+$'
options_pattern='^--\(help\|dry_run\|mode\|cluster_partition\|cluster_cpus\|cluster_memory\|cluster_exclusive\|cluster_jobs\|matlab_name\|matlab_setup\|do_not_delete_temp\)$'
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
        $@
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
            
            matlab_name)
                if [ "$#" -lt 1 ]; then
                    /bin/echo "Option --${option} requires one argument" >&2
                    exit 1
                fi
                if [ `expr length "\`expr match \"$1\" \"$integer_pattern\"\`"` -eq "0" ]; then
                    /bin/echo "Invalid matlab_name ${arg0}" >&2
                    exit 1
                fi
                matlab_name="$1"
                shift
                ;;
            
            matlab_setup)
                if [ "$#" -lt 1 ]; then
                    /bin/echo "Option --${option} requires one argument" >&2
                    exit 1
                fi
                if [ `expr length "\`expr match \"$1\" \"$integer_pattern\"\`"` -eq "0" ]; then
                    /bin/echo "Invalid matlab_setup ${arg0}" >&2
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
            
            do_not_delete_temp)
                delete_temp=0
                ;;
            
            *)
                /bin/echo "Invalid option --${option}" >&2
                exit 1
                ;;
        esac
        
    elif [ `expr length "$application_to_run"` -eq "0" ]; then
        if [ `expr match "$1" "$application_pattern"` -eq "0" ]; then

            /bin/echo "Invalid application ${arg0}" >&2
            exit 1
        fi
        application_to_run="$1"
        shift
        
    else
        /bin/echo "Too many arguments" >&2
        exit 1
    fi
done


matlab_command="${matlab_setup} ; time ${matlab_name} -nosplash -nodesktop -singleCompThread -r \"try; setup; ${application_to_run}; catch err; disp(getReport(err, 'extended')); end; exit\""
case "$mode" in
    local)
        execute_command "${matlab_command}"
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
${matlab_command} 1>"${slurm_stdout}" 2>"${slurm_stderr}"
EOF
        # cat "${tempfilename}"
        slurm_command="${slurm_command} --job-name=${application_to_run} ${tempfilename}"
        for i in $(seq 1 ${cluster_jobs}); do
            if (( i > 1 )); then
                sleep 1
            fi
            execute_command "${slurm_command}"
        done
        if [ "$delete_temp" = 1 ]; then
            rm "${tempfilename}"
        fi
        ;;
        
    *)
        /bin/echo "Invalid mode $mode" >&2
        exit 1
        ;;
esac

