#!/usr/bin/python3
'''
Generate multiple unique geometries, run simulations in them, and analyze the results for a given MCell- or Virtual Cell-format simulation description.

The VCML file should use the finite volume (deterministic PDE) solver and include biochemistry in an arbitrary geometry with compartments "EC" (extracellular region), "CP" (cytoplasm), or "NU" (nucleus).

Taraz Buck
2021-07-16

Copyright Murphy Lab 2021
'''

import sys
import traceback
import os
import glob
import re
import math
import argparse
import pathlib
import subprocess
from subprocess import run
import shlex
import time
from datetime import datetime
from datetime import timedelta
import pprint
# import numpy as np
# import xml.etree.ElementTree as ET
# import io
# import zlib
# import uuid
from itertools import chain
from vcml2fvinput import vcml2fvinput

from argparse_util import remove_prepended_arguments, existing_glob, existing_path, positive_int, bool_strict, nonnegative_int, positive_float
import biochemical_simulation_files as bsf


# cluster_mode_option = 'mode'
cluster_mode_option = 'cluster_mode'
n_images_option = 'n_images_to_synthesize'
# n_images_option = 'n_images'

def round(x):
    return int(x + 0.5)

def get_bash_path():
    # Requires Python 3.7
    # result = run('which', 'bash', capture_output=True, encoding='utf8')
    result = run(['which', 'bash'], stdout=subprocess.PIPE, encoding='utf8')
    if result.returncode != 0:
        raise FileNotFoundError('bash is required')
    bash_path = result.stdout.strip()
    return bash_path

def generate_and_simulate_info(**kw):
    '''
    Generate paths, patterns, and other values used in `generate_and_simulate` and other files that import this one.
    '''
    input_filename = str(kw['reaction_network_file'])
    input_filename_root, input_filename_ext = os.path.splitext(input_filename)
    input_filename_head, input_filename_tail = os.path.split(input_filename)
    input_filename_tail_root, input_filename_tail_ext = os.path.splitext(input_filename_tail)
    reaction_network_name = input_filename_tail_root.replace('.*', '')
    dry_run = kw['dry_run']
    downsampling = kw['downsampling']
    vcml_relative_downsampling = kw['vcml_relative_downsampling']
    synthesis = kw['synthesis']
    simulation_end_time = kw['simulation_end_time']
    #  = kw['']
    
    """
    downsampling2 = float(kw['downsampling'])
    if kw['synthesis'] == 'all':
        downsampling2 *= 1/4
    """
    
    data_set_id = bsf.generate_and_simulate_generic_data_set_id(reaction_network_name, downsampling, synthesis, simulation_end_time)
    data_set_id_vcml = bsf.generate_and_simulate_generic_data_set_id(reaction_network_name, downsampling * vcml_relative_downsampling, synthesis, simulation_end_time)
    
    mcell_simulations_directory_pattern = os.path.join(kw['output_dir'], bsf.generate_and_simulate_generic_directory_name_pattern(reaction_network_name, downsampling, synthesis, simulation_end_time))
    mcell_simulations_viz_data_directory_pattern = os.path.join(kw['output_dir'], bsf.generate_and_simulate_generic_mcell_output_pattern(reaction_network_name, downsampling, synthesis, simulation_end_time))
    vcml_simulations_directory_pattern = os.path.join(kw['output_dir'], bsf.generate_and_simulate_generic_directory_name_pattern(reaction_network_name, downsampling * vcml_relative_downsampling, synthesis, simulation_end_time))
    
    mcell_simulations_filename_main_pattern = os.path.join(mcell_simulations_directory_pattern, 'cell1', 'cell.main.mdl')
    mcell_simulations_filename_temp_pattern = mcell_simulations_directory_pattern + '.tmp'
    
    vcml_simulations_filename_main_pattern = os.path.join(vcml_simulations_directory_pattern, 'cell1', 'cell.vcml')
    vcml_simulations_filename_temp_pattern = vcml_simulations_directory_pattern + '.tmp'
    
    mcell_data_set_id = f'{data_set_id}_mcell'
    vcell_data_set_id = f'{data_set_id_vcml}_vcell'
    
    """
    print("kw['output_dir']", kw['output_dir']) # Debug
    print('input_filename_tail_root', input_filename_tail_root) # Debug
    print('reaction_network_name', reaction_network_name) # Debug
    print('downsampling2', downsampling2) # Debug
    print("kw['simulation_end_time']", kw['simulation_end_time']) # Debug
    print('simulations_directory_pattern', simulations_directory_pattern) # Debug
    print('mcell_simulations_filename_pattern', mcell_simulations_filename_pattern) # Debug
    print('mcell_simulations_filename_main_pattern', mcell_simulations_filename_main_pattern) # Debug
    raise NotImplementedError # Debug
    """
    
    def call_if_not_dry_run(func, *args, **kw2):
        kw2 = dict(kw2)
        # print('kw2'); print(kw2) # Debug
        if dry_run:
            arg_str = str(', '.join([repr(x) for x in args] + [f"'{k}': {repr(v)}" for k, v in kw2.items()]))
            # arg_str_limit = 60
            arg_str_limit = None
            if arg_str_limit is not None and len(arg_str) > arg_str_limit + 4:
                arg_str_trunc = arg_str[:arg_str_limit] + ' ...'
            else:
                arg_str_trunc = arg_str
            print(f'Dry run, would have called {func}({arg_str_trunc})')
            result = None
        else:
            result = func(*args, **kw2)
        return result
    def run_if_not_dry_run(args, capture_output=False, **kw2):
        """
        if not any([isinstance(args, x) for x in [list, tuple]):
            args = (args,)
        """
        # Requires Python 3.7
        # kw2['capture_output'] = capture_output
        if capture_output:
            kw2['stdout'] = subprocess.PIPE
            kw2['stderr'] = subprocess.PIPE
            kw2['encoding'] = 'utf8'
        return call_if_not_dry_run(run, args, **kw2)
    sbatch_output_program = re.compile(r'Submitted batch job ([0-9]+)')
    def get_sbatch_job_ids(job_submission_command_result):
        result = [sbatch_output_program.fullmatch(x) for x in job_submission_command_result.stdout.strip().splitlines()]
        result = [x for x in result if x is not None]
        return set([int(x.group(1)) for x in result])
    def get_squeue_job_ids():
        squeue_result = run_if_not_dry_run(['squeue', '-p', kw['cluster_partition']], capture_output=True)
        squeue_result_job_ids = set([int(x.strip().split()[0]) for x in squeue_result.stdout.splitlines() if not x.strip().startswith('JOBID')])
        return squeue_result_job_ids
    
    wait_time_seconds = 5
    # print_wait_time_seconds_initial = 5 * 60
    # print_wait_time_scale = math.pow(2, 1/10)
    print_wait_time_seconds_initial = 10
    print_wait_time_scale = math.pow(2, 1/4)
    def wait_for_jobs(job_ids):
        if kw[cluster_mode_option] == 'slurm':
            slurm_job_ids = set(job_ids)
            if dry_run:
                print(f'Dry run, would have waited for {len(slurm_job_ids): 4d} cluster jobs to finish.')
            else:
                # Wait for all jobs to finish
                print_wait_time_seconds = float(print_wait_time_seconds_initial)
                start_time = datetime.utcnow()
                print_last_notification_time = start_time
                time_passed = timedelta(seconds=0)
                time_passed_seconds = time_passed.total_seconds()
                while True:
                    # Check if slurm jobs are still running
                    squeue_result_job_ids = set(get_squeue_job_ids())
                    squeue_result_job_ids_ours = slurm_job_ids & squeue_result_job_ids
                    if len(squeue_result_job_ids_ours) == 0:
                        print(f'Finished waiting for {len(slurm_job_ids): 4d} cluster jobs to finish. {str(time_passed)[:-6]} total.')
                        break
                    time.sleep(wait_time_seconds)
                    current_time = datetime.utcnow()
                    time_passed = current_time - start_time
                    time_passed_seconds = time_passed.total_seconds()
                    time_passed_since_last_notification = current_time - print_last_notification_time
                    time_passed_since_last_notification_seconds = time_passed_since_last_notification.total_seconds()
                    if round(time_passed_since_last_notification_seconds) >= print_wait_time_seconds:
                        print_last_notification_time = current_time
                        print_wait_time_seconds *= print_wait_time_scale
                        # print(f'Waiting for {len(squeue_result_job_ids_ours): 4d} out of {len(slurm_job_ids): 4d} cluster jobs to finish. {str(time_passed)[:-6]} so far, next notification in {str(timedelta(seconds=print_wait_time_seconds))[:-6]}.')
                        print(f'Waiting for {len(squeue_result_job_ids_ours): 4d} out of {len(slurm_job_ids): 4d} cluster jobs to finish. {str(time_passed)[:-6]} so far, next notification in {str(timedelta(seconds=print_wait_time_seconds))[:-6]}.')
                        # print(f'Waiting for {len(squeue_result_job_ids_ours): 4d} out of {len(slurm_job_ids): 4d} cluster jobs to finish. {str(time_passed)[:-6]} so far, next notification in {print_wait_time_seconds:4f} seconds.')
        else:
            raise NotImplementedError
    
    output_vcml = input_filename_ext.lower() == '.vcml'
    output_mcell = input_filename_ext.lower() in ['.mdl', '.mcell']
    
    info = {}
    info['input_filename'] = input_filename
    info['input_filename_root'] = input_filename_root
    info['input_filename_ext'] = input_filename_ext
    info['input_filename_head'] = input_filename_head
    info['input_filename_tail'] = input_filename_tail
    info['input_filename_tail_root'] = input_filename_tail_root
    info['input_filename_tail_ext'] = input_filename_tail_ext
    info['reaction_network_name'] = reaction_network_name
    info['downsampling'] = downsampling
    info['vcml_relative_downsampling'] = vcml_relative_downsampling
    info['synthesis'] = synthesis
    info['simulation_end_time'] = simulation_end_time
    info['dry_run'] = dry_run
    info['data_set_id'] = data_set_id
    info['data_set_id_vcml'] = data_set_id_vcml
    info['mcell_simulations_directory_pattern'] = mcell_simulations_directory_pattern
    info['mcell_simulations_viz_data_directory_pattern'] = mcell_simulations_viz_data_directory_pattern
    info['vcml_simulations_directory_pattern'] = vcml_simulations_directory_pattern
    info['mcell_simulations_filename_main_pattern'] = mcell_simulations_filename_main_pattern
    info['mcell_simulations_filename_temp_pattern'] = mcell_simulations_filename_temp_pattern
    info['vcml_simulations_filename_main_pattern'] = vcml_simulations_filename_main_pattern
    info['vcml_simulations_filename_temp_pattern'] = vcml_simulations_filename_temp_pattern
    info['mcell_data_set_id'] = mcell_data_set_id
    info['vcell_data_set_id'] = vcell_data_set_id
    info['call_if_not_dry_run'] = call_if_not_dry_run
    info['run_if_not_dry_run'] = run_if_not_dry_run
    # info['sbatch_output_program'] = sbatch_output_program
    info['get_sbatch_job_ids'] = get_sbatch_job_ids
    # info['get_squeue_job_ids'] = get_squeue_job_ids
    info['wait_for_jobs'] = wait_for_jobs
    info['output_vcml'] = output_vcml
    info['output_mcell'] = output_mcell
    # info[''] = 

    print('info') # Debug
    for key, value in info.items(): # Debug
        print(f'    {key}', value) # Debug
    
    # info[''] = 
    return info
    
def generate_and_simulate(**kw):
    
    info = generate_and_simulate_info(**kw)
    input_filename = info['input_filename']
    input_filename_root = info['input_filename_root']
    input_filename_ext = info['input_filename_ext']
    input_filename_head = info['input_filename_head']
    input_filename_tail = info['input_filename_tail']
    input_filename_tail_root = info['input_filename_tail_root']
    input_filename_tail_ext = info['input_filename_tail_ext']
    reaction_network_name = info['reaction_network_name']
    downsampling = info['downsampling']
    vcml_relative_downsampling = info['vcml_relative_downsampling']
    synthesis = info['synthesis']
    simulation_end_time = info['simulation_end_time']
    dry_run = info['dry_run']
    data_set_id = info['data_set_id']
    data_set_id_vcml = info['data_set_id_vcml']
    mcell_simulations_directory_pattern = info['mcell_simulations_directory_pattern']
    mcell_simulations_viz_data_directory_pattern = info['mcell_simulations_viz_data_directory_pattern']
    vcml_simulations_directory_pattern = info['vcml_simulations_directory_pattern']
    mcell_simulations_filename_main_pattern = info['mcell_simulations_filename_main_pattern']
    mcell_simulations_filename_temp_pattern = info['mcell_simulations_filename_temp_pattern']
    vcml_simulations_filename_main_pattern = info['vcml_simulations_filename_main_pattern']
    vcml_simulations_filename_temp_pattern = info['vcml_simulations_filename_temp_pattern']
    mcell_data_set_id = info['mcell_data_set_id']
    vcell_data_set_id = info['vcell_data_set_id']
    call_if_not_dry_run = info['call_if_not_dry_run']
    run_if_not_dry_run = info['run_if_not_dry_run']
    # sbatch_output_program = info['sbatch_output_program']
    get_sbatch_job_ids = info['get_sbatch_job_ids']
    # get_squeue_job_ids = info['get_squeue_job_ids']
    wait_for_jobs = info['wait_for_jobs']
    output_vcml = info['output_vcml']
    output_mcell = info['output_mcell']
    #  = info['']
    
    # matlab_option_prefixes = {}
    # matlab_option_prefixes['']
    
    # Generate geometries and save as MCell MDL files or Virtual Cell VCML files using CellOrganizer
    # TODO Check for existence of all output and no lock files. If so, skip generation command unless recomputation requested.
    '''
    for i in $(seq 1 4) ; do /bin/bash run_application.sh --mode slurm --cluster_partition model1 --cluster_exclusive --cluster_memory 24576 --matlab_setup "module load matlab-9.7" generate_simulation_instances_min ; sleep 2 ; done
    '''
    generation_matlab_command_lines = []
    
    generation_args_kw_not_for_matlab = set()
    generation_args_kw_not_for_matlab.add('overwrite')
    generation_args_kw_not_for_matlab.add('dry_run')
    generation_args_kw_not_for_matlab.add(cluster_mode_option)
    generation_args_kw_not_for_matlab.add('cluster_partition')
    generation_args_kw_not_for_matlab.add('matlab_setup')
    generation_args_kw_not_for_matlab.add('generation_cluster_jobs')
    generation_args_kw_not_for_matlab.add('generation_cluster_memory')
    generation_args_kw_not_for_matlab.add('generation_cluster_exclusive')
    generation_args_kw_not_for_matlab.add('simulation_cluster_jobs')
    generation_args_kw_not_for_matlab.add('simulation_cluster_cpus')
    generation_args_kw_not_for_matlab.add('simulation_cluster_memory')
    generation_args_kw_not_for_matlab.add('simulation_cluster_exclusive')
    
    generation_matlab_command_lines.append('options = struct()')
    
    # Convert options for `generate_simulation_instances_generic.m` to Matlab format
    for key, val in kw.items():
        if key in generation_args_kw_not_for_matlab:
            continue
        if key == 'translations':
            val2 = zip(val[::2], val[1::2])
            val2 = '{' + '; '.join(f'{repr(x)}, {repr(y)}' for x, y in val2) + '}'
            generation_matlab_command_lines.append(f'options.{key} = {val2}')
        else:
            use_repr = True
            if isinstance(val, pathlib.Path):
                val = str(val)
            if isinstance(val, bool):
                val = str(val).lower()
                use_repr = False
            generation_matlab_command_lines.append(f'options.{key} = {repr(val) if use_repr else val}')
    
    generation_matlab_command_lines.append(f'options.output_vcml = {repr(output_vcml).lower()}')
    generation_matlab_command_lines.append(f'options.output_mcell = {repr(output_mcell).lower()}')
    # generation_matlab_command_lines.append(f'options. = {kw['']}')
    generation_matlab_command_lines.append('generate_simulation_instances_generic(options)')
    generation_matlab_command = '; '.join(generation_matlab_command_lines)
    
    """
    generation_command = []
    generation_command.extend([get_bash_path()])
    generation_command.extend(['-c'])
    """
    generation_bash_command = []
    generation_bash_command.extend([get_bash_path()])
    generation_bash_command.extend([os.path.join(kw['cellorganizer'], 'run_matlab_script.sh')])
    generation_bash_command.extend(['--mode', kw[cluster_mode_option]])
    if kw[cluster_mode_option] == 'slurm':
        generation_bash_command.extend(['--cluster_partition', kw['cluster_partition']])
        if kw['generation_cluster_exclusive']:
            generation_bash_command.extend(['--cluster_exclusive'])
        generation_bash_command.extend(['--cluster_memory', kw['generation_cluster_memory']])
        generation_bash_command.extend(['--cluster_jobs', kw['generation_cluster_jobs']])
    matlab_cd_setup = f"cd {kw['cellorganizer']}"
    if len(kw['matlab_setup']) > 0:
        matlab_cd_setup = f"{matlab_cd_setup} ; {kw['matlab_setup']}"
    generation_bash_command.extend(['--matlab_setup', matlab_cd_setup])
    generation_bash_command.extend([generation_matlab_command])
    generation_bash_command = [repr(x).lower() if isinstance(x, bool) else x for x in generation_bash_command]
    generation_bash_command = [x if isinstance(x, str) else repr(x) for x in generation_bash_command]
    generation_command = [x if isinstance(x, str) else repr(x) for x in generation_bash_command]
    
    """
    if not os.path.exists():
        os.mkdir()
    """
    mcell_main_files = glob.glob(mcell_simulations_filename_main_pattern)
    mcell_temp_files = glob.glob(mcell_simulations_filename_temp_pattern)
    should_generate_mcell = output_mcell and (len(mcell_main_files) < kw[n_images_option] or len(mcell_temp_files) > 0)
    vcml_main_files = glob.glob(vcml_simulations_filename_main_pattern)
    vcml_temp_files = glob.glob(vcml_simulations_filename_temp_pattern)
    should_generate_vcml = output_vcml and (len(vcml_main_files) < kw[n_images_option] or len(vcml_temp_files) > 0)
    if should_generate_mcell or should_generate_vcml:
        generation_command_result = run_if_not_dry_run(generation_command, capture_output=True)
        """
        if kw[cluster_mode_option] == 'slurm':
            raise NotImplementedError('Output is not captured')
        """
        if not dry_run:
            print()
            print('generation_command_result.stdout:')
            print(generation_command_result.stdout)
            print()
            print('generation_command_result.stderr:')
            print(generation_command_result.stderr)
            print()
            # raise NotImplementedError # Debug
            
            if generation_command_result.returncode != 0:
                sys.exit('generation_command failed')
                # raise subprocess.CalledProcessError('generation_command failed')
            
            if kw[cluster_mode_option] == 'slurm':
                wait_for_jobs(get_sbatch_job_ids(generation_command_result))
    
    # Run MCell simulations
    if output_mcell:
        # raise NotImplementedError # Debug
        '''
        ~/repos/cellorganizer3/applications/generate_simulation_instances_min/img_CBExMinScaled3_20min_0.50_[0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9]_4000sec
        /bin/bash ~/repos/cellorganizer3/run_mcell_simulations.sh --mode slurm --cluster_partition model1 --cluster_cpus 1 --cluster_memory 2048 --cluster_jobs 16 --seed_offset 500578 ~/repos/cellorganizer3/applications/generate_simulation_instances_min/img_CBExMinScaled3_20min_0.50_[0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9]_4000sec
        '''
        simulation_command = []
        simulation_command.extend([get_bash_path()])
        simulation_command.extend([os.path.join(kw['cellorganizer'], 'run_mcell_simulations.sh')])
        simulation_command.extend(['--mode', kw[cluster_mode_option]])
        if kw[cluster_mode_option] == 'slurm':
            simulation_command.extend(['--cluster_partition', kw['cluster_partition']])
            if kw['simulation_cluster_exclusive']:
                simulation_command.extend(['--cluster_exclusive'])
            simulation_command.extend(['--cluster_cpus', kw['simulation_cluster_cpus']])
            simulation_command.extend(['--cluster_memory', kw['simulation_cluster_memory']])
            simulation_command.extend(['--cluster_jobs', kw['simulation_cluster_jobs']])
        simulation_command.extend(['--seed_offset', kw['simulation_seed_offset']])
        # simulation_command.extend(['--overwrite_existing_results', int(kw['overwrite'])])
        if not kw['overwrite']:
            simulation_command.extend(['--keep_existing_results'])
        simulation_command.extend([mcell_simulations_directory_pattern])
        # simulation_command.extend([])
        simulation_command = [x if isinstance(x, str) else repr(x) for x in simulation_command]
        """
        if kw[cluster_mode_option] == 'slurm':
            raise NotImplementedError('Output is not captured')
        """
        # raise NotImplementedError # Debug
        print('simulation_command'); print(simulation_command) # Debug
        # raise NotImplementedError # Debug
        simulation_command_result = run_if_not_dry_run(simulation_command, capture_output=True)
        if not dry_run:
            print()
            print('simulation_command_result.stdout:')
            print(simulation_command_result.stdout)
            print()
            print('simulation_command_result.stderr:')
            print(simulation_command_result.stderr)
            print()
            
            if simulation_command_result.returncode != 0:
                sys.exit('simulation_command failed')
                # raise subprocess.CalledProcessError('simulation_command failed')
            
            # Wait
            if kw[cluster_mode_option] == 'slurm':
                wait_for_jobs(get_sbatch_job_ids(simulation_command_result))
    
    if output_vcml:
        # Convert Virtual Cell files to FiniteVolume_x64 input format using vcml2fvinput
        raise NotImplementedError
        for generated_vcml_filename in generated_vcml_filenames:
            generated_vcml_base_filename = os.path.abspath(generated_vcml_filename)
            generated_vcml_fvinput_filename = generated_vcml_base_filename + '.fvinput'
            vcml2fvinput(generated_vcml_filename, generated_vcml_fvinput_filename)
        
        # Run Virtual Cell simulations
        raise NotImplementedError
        vcell_simulations_filename_pattern = '~/repos/cellorganizer3/applications/generate_simulation_instances_min/img_CBExMinScaled3_20min_0.50_[0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9]_4000sec'
        run_if_not_dry_run([get_bash_path(), os.path.join(kw['cellorganizer'], 'run_mcell_simulations.sh'), '--mode', 'slurm', '--cluster_partition', 'model1', '--cluster_cpus', '1', '--cluster_memory', '2048', '--cluster_jobs', '16', vcell_simulations_filename_pattern])
        
        # Wait
        raise NotImplementedError
    
    """
    # Analyze results
    raise NotImplementedError
    # run()
    
    # Wait
    raise NotImplementedError
    
    # Create figures and tables
    raise NotImplementedError
    # run()
    """
    

def n_vesicles_format(x):
    x = [positive_int(y) for y in x]
    if not (len(x) == 0 or (len(x) == 2 and x[0] >= 0 and x[1] >= x[0])):
        raise argparse.ArgumentTypeError('n_vesicles must be an empty list or a list of two non-negative integers with the second entry greater than or equal to the first')
    return x


if __name__ == '__main__':
    
    # Parse options
    substitution_key_set = bsf.cellorganizer_substitution_key_set
    path_with_substitution = bsf.cellorganizer_path_with_substitution
    paths_with_substitution = bsf.cellorganizer_paths_with_substitution
    script_with_substitution = bsf.cellorganizer_path_with_substitution
    
    parser = argparse.ArgumentParser(description='Generate multiple unique geometries, run simulations in them, and analyze the results for a given MCell- or Virtual Cell-format simulation description.')
    
    parser.add_argument('reaction_network_file', type=path_with_substitution, help='Either a path of a single VCML file or a glob pattern for a single collection of MCell MDL files')
    
    parser.add_argument('output_dir', type=pathlib.Path, help='The path of a directory into which intermediate and final results will be written')
    
    parser.add_argument('--overwrite', action='store_true', help='Overwrite existing intermediate files instead of failing with error (only applies to MCell simulation output for now)')
    
    parser.add_argument('--dry_run', action='store_true', help='Print actions instead of taking them')
    
    parser.add_argument('--cellorganizer', type=existing_path, default='..', help="The path of a CellOrganizer installation")
    
    parser.add_argument(f'--{cluster_mode_option}', type=str, choices=['local', 'slurm'], default='local', help='Whether to run jobs in this process or on a Slurm cluster')
    
    # cluster_partition_default = 'pool1'
    cluster_partition_default = 'model1'
    parser.add_argument('--cluster_partition', type=str, default=cluster_partition_default, help='Cluster partition on which to queue jobs')
    
    # matlab_setup_default = 'module load matlab-9.7'
    # matlab_setup_default = '/bin/bash /etc/bashrc ; module load matlab-9.7'
    matlab_setup_default = '/bin/bash /etc/bashrc ; source "{{cellorganizer}}/module_if_available.sh" ; module_if_available load matlab-9.7'
    parser.add_argument('--matlab_setup', type=script_with_substitution, default=matlab_setup_default, help='Bash command to set up Matlab')
    
    parser.add_argument('--generation_cluster_jobs', type=positive_int, default=4, help='Number of geometry generation jobs to run')
    
    parser.add_argument('--generation_cluster_memory', type=positive_int, default=24576, help='Memory required for geometry generation job in MB')
    
    parser.add_argument('--generation_cluster_exclusive', action='store_true', help='Run only one generation job per node')
    
    # TODO Add `--simulation_cluster_nodes` argument to fill a number of nodes or a list of specified nodes
    parser.add_argument('--simulation_cluster_jobs', type=positive_int, default=16, help='Number of geometry simulation jobs to run')
    
    parser.add_argument('--simulation_cluster_cpus', type=positive_int, default=1, help='CPUs required for simulation job in MB')
    
    parser.add_argument('--simulation_cluster_memory', type=positive_int, default=2048, help='Memory required for simulation job in MB')
    
    parser.add_argument('--simulation_cluster_exclusive', action='store_true', help='Run only one simulation job per node')
    
    # framework_model_default = '{{models}}/3D/spharm/lamp2.mat'
    framework_model_default = '{{models}}/3D/spharm/lamp2.demo3D52.mat'
    parser.add_argument('--framework_model', type=path_with_substitution, default=framework_model_default, help="The path of a single CellOrganizer framework model (cell and nucleus shapes)")
    
    # parser.add_argument('--vesicle_models', type=existing_path, nargs='*', default=['{{models}}/3D/tfr.mat'], help='A list of paths of CellOrganizer vesicle models')
    parser.add_argument('--vesicle_models', type=paths_with_substitution, nargs='*', default=['{{models}}/3D/tfr.mat'], help='A list of paths of CellOrganizer vesicle models')
    
    parser.add_argument('--synthesis', type=str, choices=['framework', 'all'], default='framework', help="`'all'` to synthesize framework and vesicles or `'framework'` to ignore the vesicle model")
    
    parser.add_argument('--output_image', type=bool_strict, default=True, help='Boolean. Whether to additionally produce images of generated geometries')
    
    parser.add_argument(f'--{n_images_option}', type=positive_int, default=100, help='Number of geometries to sample from `framework_model` to produce distinct simulation inputs')
    
    parser.add_argument('--n_vesicles', type=nonnegative_int, nargs='*', default=[], help="""
    Number of vesicles to sample from each model in `vesicle_models`
        * A value of `[]` samples the number of vesicles from the model's object count distribution
        * An array of length two containing the minimum and maximum number of objects to be sampled uniformly
        """.strip())
    
    parser.add_argument('--vesicle_volume_scale', type=positive_float, default=1, help='Factor by which to scale each vesicle after synthesis')
    
    parser.add_argument('--downsampling', type=positive_float, default=1/2, help='Downsampling factor to trade disk and computational requirements for geometric precision')
    
    parser.add_argument('--vcml_relative_downsampling', type=positive_float, default=1/2, help='Factor for further downsampling for writing VCML files. VCell simulation data from fewer than 100 simulations at full resolution can fill available disk space on the VCell server, but this is less of a concern when running simulations on your own cluster.')
    
    parser.add_argument('--base_seed', type=int, default=735945, help='Seed for random number generator')
    
    parser.add_argument('--n_vesicles_seed_offset', type=int, default=34188, help='Offset used to produce seed for a second random number generator for sampling the number of vesicles')
    
    parser.add_argument('--translations', type=str, nargs='*', default=['cell', 'CP', 'nuc', 'NU', 'nucleus', 'NU', 'lamp2_mat_tfr_mat', 'EN', 'CP_EC', 'PM', 'CP_EN', 'EM', 'CP_NU', 'NM'], help='List of strings of even length where the first of each pair is to be replaced by the second in compartment names in CellOrganizer-generated VCML, MCell MDL, and SBML')
    
    parser.add_argument('--simulation_end_time', type=positive_float, default=4000, help='Time at which to end simulation in seconds. Only implemented for VCell.')
    
    parser.add_argument('--simulation_default_time_step', type=positive_float, default=1e0, help='Initial time step size in seconds for selected simulation method. Only implemented for VCell.')
    
    parser.add_argument('--simulation_max_time_step', type=positive_float, default=4, help='Maximum time step size in seconds. Only implemented for VCell.')
    
    parser.add_argument('--simulation_output_time_step', type=positive_float, default=100, help='Output time step size in seconds. Output format varies by simulation method. Only implemented for VCell.')
    
    parser.add_argument('--simulation_absolute_tolerance', type=positive_float, default=1e-8, help="VCML-specific integrator's absolute error tolerance")
    
    parser.add_argument('--simulation_relative_tolerance', type=positive_float, default=1e-8, help="VCML-specific integrator's relative error tolerance")
    
    parser.add_argument('--simulation_interaction_radius', type=positive_float, default=0.03, help='MCell-specific maximum distance for bimolecular reactions to be simulated. Not implemented.')
    
    parser.add_argument('--simulation_seed_offset', type=int, default=500578, help='Offset used to produce seed for MCell random number generator')
    
    parser.add_argument('--framework_min_clearance', type=float, default=-math.inf, help='double specifying the minimum distance in μm to impose between nucleus and cell after synthesis. -inf to disable. Currently only used for framework meshes assuming corresponding vertices by instance2MCellMDL. Default is -inf.')
    
    parser.add_argument('--framework_clearance_n_max_filter_rounds', type=positive_int, default=1, help='integer specifying the number of rounds of maximum filter to apply to the projections of cell vertices onto nucleus normals among the immediate neighbors of each vertex. 0 to disable. Currently only used for framework meshes assuming corresponding vertices by instance2MCellMDL. Default is 1.')
    
    intersecting_mesh_object_policy_default = 'ignore'
    parser.add_argument('--intersecting_mesh_object_policy', type=str, default=intersecting_mesh_object_policy_default, help='''string specifying policy for checking framework and objects for intersection and whether to remove objects or reject the synthesized cell entirely. Currently untested for values other than 'ignore'. Currently only used for framework meshes assuming corresponding vertices by instance2MCellMDL. Default is 'ignore'.''')
    
    # print('dir(parser)'); pprint.pprint(dir(parser)) # Debug
    # print('dir(parser._actions)'); pprint.pprint(dir(parser._actions)) # Debug
    
    argv = sys.argv
    argv = remove_prepended_arguments(argv)
    # print('sys.argv'); print(sys.argv) # Debug
    # print('argv'); print(argv) # Debug
    # print(''); print() # Debug
    # args = parser.parse_args(argv)
    args = parser.parse_args(argv[1:])
    args_vars = vars(args)
    
    bsf.cellorganizer_path_substitutions(parser, args_vars)
    
    # raise NotImplementedError # Debug
    
    print('args_vars'); pprint.pprint(args_vars) # Debug
    # raise NotImplementedError # Debug
    
    # args_vars['dry_run'] = True # Debug
    
    n_vesicles_format(args_vars['n_vesicles'])
    
    # generate_and_simulate(**vars(args))
    generate_and_simulate(**args_vars)
