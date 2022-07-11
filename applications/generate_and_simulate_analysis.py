'''
Python script to analyze simulation output from `generate_simulation_instances_min.m`.

Configuration options are at the top of the script.

Taraz Buck
2021-10-10
'''


from enable_interactive_debugger import sys


import os
import glob
import fnmatch
from math import *
from itertools import product
import re
import xml
import zlib
import struct
import json
import argparse
import pprint

import numpy as np
import pandas as pd
import xarray as xr
import scipy as sp
import trimesh
from sklearn.neighbors import KDTree
import statsmodels.api as sm
import statsmodels.formula.api as smf

import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import seaborn as sns



base_filename = os.path.splitext(__file__)[0]
try:
    os.mkdir(base_filename)
except OSError:
    pass


import socket
import tempfile
from datetime import datetime
import time

import pathlib
from argparse_util import remove_prepended_arguments, existing_glob, existing_path, bool_strict, positive_float

from parallel_util import chunk_start, chunk_finish, chunk_start_xarray, chunk_finish_xarray, debug_raise, debug_set_trace
import biochemical_simulation_files as bsf
from biochemical_simulation_files import should_convert_to_numpy_func, filename_sort_key, filename_sorted, contrast_func, parse_UCD_file, parse_MCell_DAT_file, plot_and_save_xarray_combinations, remove_faces, parse_shape_space_coordinates_file
from xarray_util import get_xarray_filename, save_xarray, load_xarray

from generate_and_simulate import generate_and_simulate_info



    
def generate_and_simulate_analysis(**kw):
    
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
    
    # Paths
    
    # Where to find CellOrganizer output
    simulations_base_dir = args.output_dir
    
    # Path of CellOrganizer geometry relative to `simulations_base_dir`
    # mcell_geometry_prefix = 'CBExMinScaled3EN20min_0.12'
    # mcell_geometry_prefix = 'CBExMinScaled3_20min_0.50'
    # mcell_simulation_suffix = '4000sec'
    # mcell_geometry_dir_name_pattern = os.path.join(f'img_{mcell_geometry_prefix}_[0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9]_{mcell_simulation_suffix}', 'cell1')
    # mcell_geometry_dir_name_pattern = mcell_simulations_directory_pattern
    # mcell_geometry_dir_name_pattern = os.path.join(generate_and_simulate_generic_directory_name_pattern(reaction_network_name, downsampling, synthesis, simulation_end_time), 'cell1')
    # Path of MCell simulation output relative to `mcell_geometry_dir_name_pattern`
    # mcell_simulation_data_name_pattern = os.path.join('viz_data', 'seed_[0-9]*')
    # For naming statistics and plot files
    # mcell_data_set_id = f'{mcell_geometry_prefix}_{mcell_simulation_suffix}_mcell'
    # data_set_id = generate_and_simulate_generic_data_set_id(reaction_network_name, downsampling, synthesis, simulation_end_time)
    # mcell_data_set_id = f'{data_set_id}_mcell'
    
    # debug_raise() # Debug
    
    # Path of VCell simulation output relative to `simulations_base_dir`
    # Simulation output should be volume variables exported in UCD format
    # vcell_simulation_data_name_pattern = 'img_SarmaGhosh2012ForCO_0.50_4000sec'
    # data_set_id_vcell = generate_and_simulate_generic_data_set_id(reaction_network_name, downsampling * vcml_relative_downsampling, synthesis, simulation_end_time)
    # Path of CellOrganizer geometry relative to `simulations_base_dir`
    # vcell_geometry_data_name_pattern = 'img_SarmaGhosh2012ForCO_0.50_[0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9]_4000sec'
    # For naming statistics and plot files
    # vcell_data_set_id = 'SarmaGhosh2012ForCO_0.50_vcell'
    
    
    
    # Debug options
    
    
    debug_max_simulations = None
    debug_max_times = None
    debug_max_faces = None
    debug_always_overwrite_data = False
    debug_always_overwrite_features = False
    debug_always_overwrite_analysis = False
    debug_always_overwrite_plots = False
    debug_always_load_features = False
    debug_skip_computing_features = False
    
    # debug_max_simulations = 2 # Debug
    debug_max_simulations = 5 # Debug
    # debug_max_simulations = 10 # Debug
    if debug_max_simulations is not None:
        print('\nDEBUG debug_max_simulations = {0}\n'.format(debug_max_simulations));
    
    # debug_max_times = 2 # Debug
    # debug_max_times = 5 # Debug
    debug_max_times = 10 # Debug
    if debug_max_times is not None:
        print('\nDEBUG debug_max_times = {0}\n'.format(debug_max_times));
    
    # debug_max_faces = 2 # Debug
    debug_max_faces = 5 # Debug
    if debug_max_faces is not None:
        print('\nDEBUG debug_max_faces = {0}\n'.format(debug_max_faces));
    
    debug_always_overwrite_data = True
    if debug_always_overwrite_data:
        print('\nDEBUG debug_always_overwrite_data = {0}\n'.format(debug_always_overwrite_data));
    
    debug_always_overwrite_features = True
    if debug_always_overwrite_features:
        print('\nDEBUG debug_always_overwrite_features = {0}\n'.format(debug_always_overwrite_features));
    
    debug_always_overwrite_analysis = True
    if debug_always_overwrite_analysis:
        print('\nDEBUG debug_always_overwrite_analysis = {0}\n'.format(debug_always_overwrite_analysis));
    
    debug_always_overwrite_plots = True
    if debug_always_overwrite_plots:
        print('\nDEBUG debug_always_overwrite_plots = {0}\n'.format(debug_always_overwrite_plots));
    
    debug_always_load_features = True
    if debug_always_load_features:
        print('\nDEBUG debug_always_load_features = {0}\n'.format(debug_always_load_features));
    
    # debug_skip_computing_features = True
    if debug_skip_computing_features:
        print('\nDEBUG debug_skip_computing_features = {0}\n'.format(debug_skip_computing_features));
    
    
    
    # Computation options
    
    
    # enable_mcell_analysis = False
    enable_mcell_analysis = True
    
    enable_vcell_analysis = False
    # enable_vcell_analysis = True
    
    # control_method = 'none'
    # control_method = 'qr'
    control_method = 'ols'
    
    # should_compute_regression_statistics = False
    should_compute_regression_statistics = True
    
    
    independent_features_to_control_for = []
    # independent_features_to_control_for.append('*')
    for compartment_name in ['CP', 'NU']:
        independent_features_to_control_for.append(f'Volume {compartment_name}')
        independent_features_to_control_for.append(f'Surface area {compartment_name}')
        for dim_str in 'XYZ':
            independent_features_to_control_for.append(f'Center of mass {compartment_name} {dim_str}')
            independent_features_to_control_for.append(f'Size {compartment_name} {dim_str}')
    # raise NotImplementedError('Control for COG and dimensions')
    
    # independent_features_control_max_exponent = 1
    independent_features_control_max_exponent = 4
    if debug_max_simulations is not None:
        independent_features_control_max_exponent = min(independent_features_control_max_exponent, 2)
    
    space_space_dimension_names_func = lambda x: ['Shape space dimension {0:d}'.format(i) for i in range(x)]
    independent_features_with_measurable_mean = []
    independent_features_with_measurable_mean.extend(space_space_dimension_names_func(8))
    # independent_features_with_measurable_mean.append('')
    
    # relative_coefficient_method = 'const'
    relative_coefficient_method = 'std'
    
    geometric_feature_prefix = 'Geometric feature'
    predictor_feature_generic_name = 'Predictor feature'
    
    
    
    # Plot options (some untested)
    
    
    should_plot_simulations = False
    # should_plot_simulations = True
    # if not should_plot_simulations:
        # print('\nDEBUG should_plot_simulations = {0}\n'.format(should_plot_simulations));
    
    should_plot_pairs = False
    # should_plot_pairs = True
    # if not should_plot_pairs:
        # print('\nDEBUG should_plot_pairs = {0}\n'.format(should_plot_pairs));
    
    should_plot_pair_grids = False
    # should_plot_pair_grids = True
    # if not should_plot_pair_grids:
        # print('\nDEBUG should_plot_pair_grids = {0}\n'.format(should_plot_pair_grids));
    
    should_plot_pair_mean_concentrations = False
    # should_plot_pair_mean_concentrations = True
    # if not should_plot_pair_mean_concentrations:
        # print('\nDEBUG should_plot_pair_mean_concentrations = {0}\n'.format(should_plot_pair_mean_concentrations));
    
    should_plot_simulation_mean_concentrations = False
    # should_plot_simulation_mean_concentrations = True
    # if not should_plot_simulation_mean_concentrations:
        # print('\nDEBUG should_plot_simulation_mean_concentrations = {0}\n'.format(should_plot_simulation_mean_concentrations));
    
    should_plot_variable_mean_concentrations = False
    # should_plot_variable_mean_concentrations = True
    # if not should_plot_variable_mean_concentrations:
        # print('\nDEBUG should_plot_variable_mean_concentrations = {0}\n'.format(should_plot_variable_mean_concentrations));
    
    
    
    
    # Functions
    
    
    float_pattern = r'(?:[\-\+]?(?:(?:(?:(?<![0-9.])[0-9][0-9]*(?![0-9.])|[0-9]+\.[0-9]*(?![0-9])|(?<![0-9])[0-9]*\.[0-9]+)(?:[eE][\-\+]?[0-9]+)?)|[iI][nN][fF]|[nN][aA][nN]))'
    float_pattern_program = re.compile(float_pattern)
    # Test:
    float_pattern_test_input = '5 5. .5 0.5 15 15. .15 15.15 -5 -5. -.5 -0.5 -15 -15. -.15 -15.15 5e3 5.e3 .5e3 0.5e3 15e3 15.e3 .15e3 15.15e3 +5e-3 +5.e-3 +.5e-3 +0.5e-3 +15e-3 +15.e-3 +.15e-3 +15.15e-3 inf Inf INF nan NaN NAN'
    float_pattern_test_output = float_pattern_program.findall(float_pattern_test_input)
    float_pattern_program_matches = lambda x: float_pattern_program.fullmatch(x) is not None
    assert(float_pattern_test_output == float_pattern_test_input.split())
    temp1 = [float(x) for x in float_pattern_test_output]
    temp2 = [float(x) for x in float_pattern_test_input.split()]
    temp3 = [x == y or (isnan(x) and isnan(y)) for x, y in zip(temp1, temp2)]
    assert(temp3)
    
    def parse_cell_info_file(given_file_path):
        '''
        Parse output from CellOrganizer with `options.output.writeMCellMDL = true`.
        
        Assumes there is one object per mesh.
        '''
        
        filename = os.path.split(given_file_path)[1]
        geometry_name = given_file_path.split(os.path.sep)[-3]
        
        with open(given_file_path, 'r') as given_file:
            given_file_text = given_file.read()
        
        json_data = json.loads(given_file_text)
        json_data_info = json_data['geometry_info']
        
        """
        values = {}
        for model_meshes_data in filter(lambda x: x['name'] == 'modelMesh', json_data['meshData']):
            if isinstance(model_meshes_data['list'], dict):
                model_meshes_data['list'] = [model_meshes_data['list']]
            for model_mesh_data in model_meshes_data['list']:
                model_name = model_mesh_data['name']
                if model_name not in values:
                    values[model_name] = 0
                values[model_name] += 1
        
        cell_info_names = [f'Object count {model_name}' for model_name in sorted(values.keys())]
        cell_info_values = [values[model_name] for model_name in sorted(values.keys())]
        
        result = pd.Series(cell_info_values, name=geometry_name, index=pd.Index(cell_info_names, name=geometric_feature_prefix), dtype='uint32')
        """
        
        models_meshes_data = list(filter(lambda x: x['name'] == 'modelMesh', json_data['meshData']))
        models_names = sorted([x['name'] for x in models_meshes_data])
        # compartment_names = json_data_info['all_compartment_names']
        compartment_names = json_data_info['names']
        # print('Excluding EC')
        compartment_names = [x for x in compartment_names if x != 'EC']
        # Surface areas and volumes possibly incorrect
        # compartment_volumes = json_data_info['all_compartment_volumes']
        # compartment_exclusive_volumes = json_data_info['all_compartment_exclusive_volumes']
        # compartment_surface_areas = json_data_info['all_compartment_surface_areas']
        
        # Surface areas and volumes possibly incorrect, load mesh and compute geometric features
        compartment_volumes = [np.nan for i in range(len(compartment_names))]
        compartment_surface_areas = [np.nan for i in range(len(compartment_names))]
        compartment_surface_coms = [(np.nan,) * 3 for i in range(len(compartment_names))]
        compartment_surface_sizes = [(np.nan,) * 3 for i in range(len(compartment_names))]
        # compartment_coms = [[None] * 3 for i in range(len(compartment_names))]
        # raise NotImplementedError('Matlab generation of the JSON file is WIP and untested')
        for i, compartment_name in enumerate(compartment_names):
            # mesh_filename = os.path.join(os.path.split(given_file_path)[0], f'{compartment_name}.obj')
            mesh_filename = os.path.join(os.path.split(given_file_path)[0], f'mcell_{compartment_name}.obj')
            if os.path.exists(mesh_filename):
                mesh = trimesh.load_mesh(mesh_filename)
                compartment_surface_areas[i] = mesh.area
                compartment_volumes[i] = mesh.volume
                if mesh.is_watertight:
                    compartment_surface_coms[i] = mesh.center_mass
                compartment_surface_sizes[i] = mesh.bounds[1]
                # compartment_coms[i] = mesh.center_mass
        compartment_surface_coms = np.array(compartment_surface_coms)
        compartment_surface_sizes = np.array(compartment_surface_sizes)
        
        compartment_name_to_index_dict = {y: x for x, y in enumerate(compartment_names)}
        
        values = {}
        
        models_features_data = []
        models_features_data.append((compartment_names, compartment_volumes, 'Volume'))
        models_features_data.append((compartment_names, compartment_surface_areas, 'Surface area'))
        models_features_data.append((compartment_names, compartment_surface_coms, 'Center of mass'))
        models_features_data.append((compartment_names, compartment_surface_sizes, 'Size'))
        # models_features_data.append((compartment_names, compartment_exclusive_volumes, 'Exclusive volume'))
        models_features_data.append((compartment_names, [np.pi**(1/3) * (6 * x)**(2/3) / y for x, y in zip(compartment_volumes, compartment_surface_areas)], 'Sphericity'))
        temp = {x: compartment_volumes[compartment_name_to_index_dict[x]] for x in compartment_names}
        if 'EN' not in temp.keys():
            temp['EN'] = np.nan
        models_features_data.append((['NU CP'], [temp['NU'] / temp['CP']], 'Volume ratio'))
        models_features_data.append((['EN CP-NU'], [temp['EN'] / (temp['CP'] - temp['NU'])], 'Volume ratio'))
        
        models_features_data = [(x, np.array(y), z) for x, y, z in models_features_data]
        
        for object_names, feature_values, feature_prefix in models_features_data:
            if feature_values.ndim == 1:
                for object_name, feature_value in zip(object_names, feature_values):
                    feature_name = f'{feature_prefix} {object_name}'
                    values[feature_name] = feature_value
            else:
                if feature_values.shape[1] != 3:
                    raise NotImplementedError
                for dim, dim_str in enumerate('XYZ'):
                    for object_name, feature_value in zip(object_names, feature_values[:, dim]):
                        feature_name = f'{feature_prefix} {object_name} {dim_str}'
                        values[feature_name] = feature_value
        
        for model_meshes_data in models_meshes_data:
            if isinstance(model_meshes_data['list'], dict):
                model_meshes_data['list'] = [model_meshes_data['list']]
            for model_mesh_data in model_meshes_data['list']:
                model_name = model_mesh_data['name']
                object_count_name = f'Object count {model_name}'
                if object_count_name not in values:
                    values[object_count_name] = 0
                values[object_count_name] += 1
        
        values_names = sorted(values.keys())
        values = [[values[x] for x in values_names]]
        result = xr.DataArray(values, dims=['Geometry', geometric_feature_prefix], coords={'Geometry': [geometry_name], geometric_feature_prefix: values_names})
        
        # cell_info_names = [f'Object count {model_name}' for model_name in sorted(values.keys())]
        # cell_info_values = [values[model_name] for model_name in sorted(values.keys())]
        
        # result = pd.Series(cell_info_values, name=geometry_name, index=pd.Index(cell_info_names, name=), dtype='uint32')
        
        return result
    
    
    def load_geometric_coordinates(coords_base_dir, coords_dir_pattern=None):
        coords_base_dir = os.path.expanduser(coords_base_dir)
        if coords_dir_pattern is not None:
            coords_base_dir = os.path.join(coords_base_dir, coords_dir_pattern)
        coords_dir_list = glob.glob(os.path.join(coords_base_dir))
        coords_dir_list = list(coords_dir_list)
        coords_dir_list = map(lambda x: os.path.join(x, 'cell1') if not x.endswith('cell1') else x, coords_dir_list)
        coords_dir_list = filename_sorted(coords_dir_list)
        coords_dir_list = list(coords_dir_list)
    
        all_shape_space_coordinates = []
        all_cell_info = []
    
        for coords_dir in coords_dir_list:
            # Parse files
            
            coords_path = os.path.join(coords_dir, 'shape_space_coords.txt');
            cell_info_path = os.path.join(coords_dir, 'cell.info.json');
            
            if not all([os.path.exists(x) for x in [coords_path, cell_info_path]]):
                continue
            
            shape_space_coordinates = parse_shape_space_coordinates_file(coords_path, index_name=geometric_feature_prefix)
            cell_info = parse_cell_info_file(cell_info_path)
            
            if shape_space_coordinates is None or cell_info is None:
                continue
            
            all_shape_space_coordinates.append(shape_space_coordinates)
            all_cell_info.append(cell_info)
            
        all_shape_space_coordinates = pd.concat(all_shape_space_coordinates, axis=1).T
        all_shape_space_coordinates.index.name = 'Geometry'
        all_shape_space_coordinates.columns.name = geometric_feature_prefix
        
        # all_cell_info = pd.concat(all_cell_info, axis=1).T
        # all_cell_info.index.name = 'Geometry'
        # all_cell_info.columns.name = geometric_feature_prefix
        all_cell_info = xr.concat(all_cell_info, dim='Geometry')
        # raise NotImplementedError
        
        geometric_coordinates = []
        geometric_coordinates.append(all_shape_space_coordinates)
        geometric_coordinates.append(all_cell_info)
        
        geometric_coordinates = xr.concat([xr.DataArray(x) for x in geometric_coordinates], dim=geometric_feature_prefix)
        geometric_coordinates = geometric_coordinates.rename({'Geometry': 'Geometry name'})
        
        return geometric_coordinates
    
    
    
    def ucd_filename_to_simulation_id(given_filename):
        return filename_sort_key(given_filename[1])
    
    
    def compute_UCD_preprocessing(given_all_UCD_info, given_save_mesh=False, given_ignore_constant_faces=False):
        '''
        '''
        
        given_UCD_info_keys = filename_sorted(given_all_UCD_info.keys())
        n_given_UCD_info_keys = len(given_UCD_info_keys)
        
        
        # Process mesh
        
        given_UCD_info_key = given_UCD_info_keys[0]
        given_UCD_info = given_all_UCD_info[given_UCD_info_key]
        file_path = given_UCD_info['file_path']
        filename = given_UCD_info['filename']
        vertices = given_UCD_info['vertices']
        faces = given_UCD_info['faces']
        face_data = given_UCD_info['face_data']
        faces_types = given_UCD_info['faces_types']
        face_attributes = {'face_data': face_data}
    
        # Find faces where the concentration changes
        
        variables_faces_masks = (face_data * 0).astype('bool')
        
        for given_UCD_info_key1, given_UCD_info_key2 in zip(given_UCD_info_keys[:-1], given_UCD_info_keys[1:]):
            face_data1 = given_all_UCD_info[given_UCD_info_key1]['face_data']
            face_data2 = given_all_UCD_info[given_UCD_info_key2]['face_data']
            if face_data1.shape != face_data2.shape:
                raise ValueError
            variables_faces_masks |= (face_data1 != face_data2)
        
        faces_mask = variables_faces_masks.any(dim='Variable')
        
        # Compute geometric measures
        # Assumes faces are voxels
        variables_faces_masks_volumes = variables_faces_masks.sum()
        given_simulation_id = ucd_filename_to_simulation_id(given_UCD_info_key)
        
        given_ignore_constant_faces2 = bool(given_ignore_constant_faces)
        if debug_max_faces is not None:
            # Debug
            # Keep the first debug_max_faces faces with changing values
            faces_mask2 = faces_mask[faces_mask]
            faces_mask2[debug_max_faces:] = False
            faces_mask[faces_mask] = faces_mask2
            given_ignore_constant_faces2 = True
        
        if given_ignore_constant_faces2:
            # Remove faces where the concentration does not change (not participating in reaction or diffusion)
            vertices, faces2, face_attributes, vertices_mask = remove_faces(vertices, faces, face_attributes, np.array(faces_mask))
            faces = faces[faces_mask, ...]
        
        if (faces_types == 'quad').all():
            # Triangulate square
            faces = np.vstack((faces[:, [0, 1, 2]], faces[:, [2, 3, 0]]))
            for key, value in face_attributes.items():
                face_attributes[key] = np.vstack((value,) * 2)
        elif (faces_types == 'hex').all():
            cube_mode = 'centroid'
            # cube_mode = 'triangulate'
            if cube_mode == 'centroid':
                # Convert cube to centroid
                vertices2 = vertices.loc[faces[:, 0], :].copy()
                for i in range(1, faces.shape[1]):
                    vertices2 += vertices.loc[faces[:, i], :]
                vertices2 /= faces.shape[1]
                vertices2 = vertices2.reset_coords('Vertex index', drop=True)
                vertices = vertices2
                faces = xr.DataArray(np.arange(vertices.shape[0])[:, None], dims=('Face index', 'Face vertex index'), coords={'Face index': np.arange(vertices.shape[0]), 'Face vertex index': [0]})
            elif cube_mode == 'triangulate':
                # Triangulate cube
                faces = np.vstack([
                    faces[:, [0, 2, 1]],
                    faces[:, [2, 0, 3]],
                    faces[:, [0, 1, 5]],
                    faces[:, [0, 5, 4]],
                    faces[:, [0, 4, 7]],
                    faces[:, [0, 7, 3]],
                    faces[:, [1, 2, 6]],
                    faces[:, [1, 6, 5]],
                    faces[:, [2, 3, 7]],
                    faces[:, [2, 7, 6]],
                    faces[:, [4, 5, 6]],
                    faces[:, [4, 6, 7]],
                    ])
                for key, value in face_attributes.items():
                    face_attributes[key] = np.vstack((value,) * 12)
        
        # Create Trimesh and save
        mesh = trimesh.Trimesh(vertices, faces, face_attributes=face_attributes, process=False)
        if given_save_mesh:
            export_file_path = os.path.join(base_filename, filename + '.obj')
            mesh.export(export_file_path)
        
        # Collect face data from all frames
        all_face_data = []
        for i, given_UCD_info_key in enumerate(given_UCD_info_keys):
            given_UCD_info = given_all_UCD_info[given_UCD_info_key]
            
            face_data = given_UCD_info['face_data']
            
            if given_ignore_constant_faces2:
                face_data = face_data.loc[faces_mask, :]
            all_face_data.append(face_data)
        
        all_face_data2 = xr.concat(all_face_data, dim=pd.Index(np.arange(n_given_UCD_info_keys), name='Time index'))
        if all_face_data2.size == 0:
            raise ValueError
        all_face_data = all_face_data2
        
        return_value = {}
        return_value['all_face_data'] = all_face_data
        return_value['vertices'] = vertices
        return_value['faces'] = faces
        return_value['faces_mask'] = faces_mask
        return return_value
    
    
    def compute_UCD_simulation_features_20200203(given_simulation_concentrations_masked):
        '''
        Returns a DataArray where the first dimension is geometry, the last dimension is feature, and other dimensions will be iterated over to build individual models.
        '''
        
        given_concentrations = given_simulation_concentrations_masked
        
        
        # Compute frame features
    
        result_features = {}
        result_features['t00_slope_mean'] = given_concentrations.differentiate('Time index', edge_order=2)[{'Time index': 0}].mean('Face index')
        final_time_prefix = 't{0:2d}'.format(given_concentrations.coords['Time index'].values[-1])
        final_minus_half_period_time_prefix = 't{0:2d}'.format(given_concentrations.coords['Time index'].values[-7])
        result_features[final_time_prefix + '_mean'] = given_concentrations[{'Time index': -1}].mean('Face index')
        result_features[final_time_prefix + '_std'] = given_concentrations[{'Time index': -1}].std('Face index')
        feature3_name = final_time_prefix + '_minus_' + final_minus_half_period_time_prefix + '_mean'
        feature3 = given_concentrations[{'Time index': -1}].reset_coords(drop=True) - given_concentrations[{'Time index': -7}].reset_coords(drop=True)
        feature3 = feature3.mean('Face index')
        result_features[feature3_name] = feature3
        
        for result_feature_name, result_feature in result_features.items():
            result_features[result_feature_name] = result_feature.reset_coords(drop=True)
        
        result_features_names = sorted(result_features.keys())
        result_features = xr.concat([result_features[x] for x in result_features_names], dim=pd.Index(result_features_names, name='Feature'))
        
        return result_features
    
    
    def mcell_seed_name_seed_func(seed_name):
        return(int(re.fullmatch(r'seed_([0-9]+)', seed_name).group(1)))
    
    def mcell_geometry_name_abbreviation_func(geometry_name):
        geometry_name = 'img_{}'.format(bsf.generate_and_simulate_generic_directory_name_geometry_seed(geometry_name))
        return geometry_name
    
    distribution_feature_funcs = {}
    distribution_feature_funcs['mean'] = lambda x, y: x.mean(dim=y)
    distribution_feature_funcs['std'] = lambda x, y: x.std(dim=y)
    distribution_feature_funcs['skew'] = lambda x, y: x.reduce(sp.stats.skew, dim=y, nan_policy='omit')
    distribution_feature_funcs['kurtosis'] = lambda x, y: x.reduce(sp.stats.kurtosis, dim=y, nan_policy='omit')
    # distribution_feature_funcs['kurtosis'] = lambda x, y: x.reduce(sp.stats.kurtosis, dim=y) if y is not None and len(y) > 0 else sp.stats.kurtosis(x)
    distribution_feature_contrast_quantile = 0.05
    distribution_feature_contrast_feature_name = 'contrast{0:0.3f}'.format(distribution_feature_contrast_quantile)
    distribution_feature_funcs[distribution_feature_contrast_feature_name] = lambda x, y: contrast_func(x, distribution_feature_contrast_quantile, y)
    
    distribution_feature_names = sorted(distribution_feature_funcs.keys())
    
    def distribution_features_func(given_features, given_distribution_feature_names, given_summary_dimension_name='Iteration', given_distribution_feature_dimension_name='Iteration feature'):
        '''
        Calculate `given_distribution_feature_names` across dimension `given_summary_dimension_name`
        '''
        result_dims = [given_distribution_feature_dimension_name]
        # result_coords_keys = [given_distribution_feature_dimension_name]
        result_coords_values = [given_distribution_feature_names]
        result_shape = tuple([len(x) for x in result_coords_values])
        if should_convert_to_numpy_func(given_features):
            given_features = xr.DataArray(given_features)
        if given_summary_dimension_name is not None:
            # result_dims = [given_distribution_feature_dimension_name] + result_dims
            given_remaining_dims = list(filter(lambda x: x != given_summary_dimension_name, given_features.dims))
            result_dims += given_remaining_dims
            # result_coords_keys = list(given_features.coords.keys()) + result_coords_keys
            # result_coords_values = list(given_features.coords.values()) + result_coords_values
            result_coords_values += [given_features.coords[x].values for x in given_remaining_dims]
            # result_shape = result_shape + tuple([len(result_coords_values) for x in result_coords_values])
            result_shape = [len(x) for x in result_coords_values]
        # result_coords = dict(zip(result_coords_keys, result_coords_values))
        result_coords = dict(zip(result_dims, result_coords_values))
        result_features = np.zeros(result_shape) + np.nan
        result_features = xr.DataArray(result_features, dims=result_dims, coords=result_coords)
        for given_distribution_feature_name in given_distribution_feature_names:
            temp = distribution_feature_funcs[given_distribution_feature_name](given_features, given_summary_dimension_name)
            result_features.loc[{given_distribution_feature_dimension_name: given_distribution_feature_name}] = temp
            if float(temp.isnull().mean()) >= 0.5:
                debug_raise()
                # debug_set_trace()
        return result_features
    
    def mcell_analysis(info, data_set_id):
        '''
        Preprocess simulation data for plotting and regression analysis
        '''
        # simulation_data_dir_list = glob.glob(info['mcell_simulations_directory_pattern'])
        simulation_data_dir_list = glob.glob(info['mcell_simulations_viz_data_directory_pattern'])
        simulation_data_dir_list = filename_sorted(simulation_data_dir_list)
        
        data_set_base_filename = os.path.join(base_filename, data_set_id)
        if debug_max_simulations is not None:
            data_set_base_filename += '_maxsim{0:06d}'.format(debug_max_simulations)
        if debug_max_times is not None:
            data_set_base_filename += '_maxt{0:06d}'.format(debug_max_times)
        if debug_max_faces is not None:
            data_set_base_filename += '_maxf{0:06d}'.format(debug_max_faces)
        
        print('len(simulation_data_dir_list) "{}"'.format(len(simulation_data_dir_list)))
        if debug_max_simulations is not None:
            simulation_data_dir_list = simulation_data_dir_list[:debug_max_simulations]
            print('len(simulation_data_dir_list) "{}"'.format(len(simulation_data_dir_list)))
        
        
        geometry_names = []
        geometry_seeds = []
        seeds = []
        
        simulations_iterations_features = []
        simulations_features = []
        simulations_features_all_included = True
        
        simulation_base_filenames = []
        
        for simulation_data_dir in simulation_data_dir_list:
            simulation_data_dir_names = simulation_data_dir.split(os.path.sep)
            simulation_geometry_full_name = simulation_data_dir_names[-4]
            simulation_geometry_seed = bsf.generate_and_simulate_generic_directory_name_geometry_seed(simulation_data_dir_names[-4])
            simulation_geometry_name = mcell_geometry_name_abbreviation_func(simulation_geometry_full_name)
            simulation_seed_name = simulation_data_dir_names[-1]
            simulation_seed = int(mcell_seed_name_seed_func(simulation_seed_name))
            simulation_base_filename = '{} {} {}'.format(data_set_base_filename, simulation_geometry_name, simulation_seed_name)
            geometry_names.append(simulation_geometry_name)
            geometry_seeds.append(simulation_geometry_seed)
            seeds.append(simulation_seed)
            simulation_base_filenames.append(simulation_base_filename)
            
            simulation_iterations_features_filename = get_xarray_filename(simulation_base_filename, 'iterations_features')
            simulations_features_filename = get_xarray_filename(simulation_base_filename, 'features')
            
            # Load data, compute features per time point
            
            can_start, final_exists, temp_name = chunk_start_xarray(simulation_base_filename, 'iterations_features')
            if debug_always_overwrite_features and final_exists:
                can_start = True
            
            if not debug_always_overwrite_features and final_exists:
                
                simulation_iterations_features = load_xarray(simulation_base_filename, ['iterations_features'])['iterations_features']
                
            elif can_start:
                
                file_path_list = glob.glob(os.path.join(simulation_data_dir, '*.dat'))
                file_path_list = filename_sorted(file_path_list)
                
                if debug_max_times is not None:
                    file_path_list = file_path_list[:debug_max_times]
                
                n_files = len(file_path_list)
                is_first_file = True
                simulations_data = {}
                # variable_names = None
                variable_names = set()
                simulations_iterations = []
                simulations_molecules_counts = []
                for j, file_path in enumerate(file_path_list):
                    
                    simulation_data = parse_MCell_DAT_file(file_path)
                    
                    filename = simulation_data['filename']
                    iteration = simulation_data['iteration']
                    
                    simulations_data[iteration] = simulation_data
                    simulations_iterations.append(iteration)
                    
                    """
                    if variable_names is None:
                    """
                    variable_names |= set(simulation_data['variable_names'])
                n_variable_names = len(variable_names)
                variable_names = sorted(variable_names)
                    
                simulations_molecules_counts = xr.DataArray(np.zeros((n_files, n_variable_names)) + np.nan, dims=('Iteration', 'Variable'), coords={'Iteration': simulations_iterations, 'Variable': variable_names}, name='molecule_count')
                
                neighbor_k = 5
                simulations_neighbor_distance_summary_features = xr.DataArray(np.zeros((n_files, n_variable_names, len(distribution_feature_names))) + np.nan, dims=('Iteration', 'Variable', 'Iteration feature'), coords={'Iteration': simulations_iterations, 'Variable': variable_names, 'Iteration feature': distribution_feature_names}, name='neighbor_distance_summary_features')
                
                simulations_molecule_coords_mean = xr.DataArray(np.zeros((n_files, n_variable_names, 3)) + np.nan, dims=('Iteration', 'Variable', 'Dimension'), coords={'Iteration': simulations_iterations, 'Variable': variable_names, 'Dimension': ['x', 'y', 'z']}, name='molecule_coords_mean')
                
                for iteration, simulation_data in simulations_data.items():
                    simulation_variable_names = simulation_data['variable_names']
                    
                    simulations_molecules_counts.loc[{'Iteration': iteration, 'Variable': simulation_variable_names}] = simulation_data['n_molecules']
                    
                    for variable_name in simulation_variable_names:
                        iteration_xyz = simulation_data['DataArray'][variable_name].sel(coord=['x', 'y', 'z'])
                        iteration_kdtree = KDTree(iteration_xyz)
                        
                        iteration_neighbor_distances, iteration_neighbor_indices = iteration_kdtree.query(iteration_xyz, k=min(neighbor_k + 1, iteration_xyz.shape[0]))
                        simulations_neighbor_distance_summary_features.loc[{'Iteration': iteration, 'Variable': variable_name}] = distribution_features_func(iteration_neighbor_distances[:, -1], distribution_feature_names, None)
                        
                        simulations_molecule_coords_mean.loc[{'Iteration': iteration, 'Variable': variable_name}] = iteration_xyz.mean(dim='index')
                
                simulation_iterations_features = {}
                simulation_iterations_features['molecule_count'] = simulations_molecules_counts
                for feature_name in simulations_neighbor_distance_summary_features.coords['Iteration feature'].values:
                    simulation_iterations_features[f'molecule_neighbor{neighbor_k}_distance_{feature_name}'] = simulations_neighbor_distance_summary_features.loc[{'Iteration feature': feature_name}].reset_coords(drop=True)
                for dim in simulations_molecule_coords_mean.coords['Dimension'].values:
                    simulation_iterations_features[f'molecule_{dim}_mean'] = simulations_molecule_coords_mean.loc[{'Dimension': dim}].reset_coords(drop=True)
                
                # raise NotImplementedError('Remove Dimension!')
                
                simulation_iterations_features = xr.concat(list(simulation_iterations_features.values()), dim=pd.Index(list(simulation_iterations_features.keys()), name='Iteration feature'))
                
                save_xarray(simulation_base_filename, iterations_features=simulation_iterations_features)
                
                chunk_finish_xarray(simulation_base_filename, 'iterations_features')
                final_exists = True
            
            if final_exists:
                # Fix nan concentration when molecule is not present
                simulation_iterations_features.loc[{'Iteration feature': 'molecule_count', 'Iteration': 0, 'Variable': simulation_iterations_features.loc[{'Iteration feature': 'molecule_count', 'Iteration': 0}].isnull()}] = 0 
                simulations_iterations_features.append(simulation_iterations_features)
            
            
            # Compute features across time
            
            can_start, final_exists, temp_name = chunk_start_xarray(simulation_base_filename, 'features')
            if debug_always_overwrite_features and final_exists:
                can_start = True
            
            if not debug_always_overwrite_features and final_exists:
                
                simulation_features = load_xarray(simulation_base_filename, ['features'])['features']
                
            elif can_start:
                
                if not debug_skip_computing_features:
                    simulation_features = []
                    for iteration_feature_name in simulation_iterations_features.coords['Iteration feature'].values:
                        simulation_feature_names = [' '.join([iteration_feature_name, x]) for x in distribution_feature_names]
                        iteration_feature = simulation_iterations_features.sel({'Iteration feature': iteration_feature_name}).reset_coords(drop=True)
                        simulation_features.append(distribution_features_func(iteration_feature, distribution_feature_names, 'Iteration'))
                    
                    # debug_raise()
                    simulation_features = xr.concat(simulation_features, dim=pd.Index(simulation_iterations_features.coords['Iteration feature'].values, name='Feature'))
                    
                    simulation_features.name = simulation_geometry_name
                
                save_xarray(simulation_base_filename, features=simulation_features)
                
                chunk_finish_xarray(simulation_base_filename, 'features')
                final_exists = True
            
            if final_exists:
                simulations_features.append(simulation_features)
            else:
                simulations_features_all_included = False
            
            # Optionally load other data
            
            if os.path.exists(simulations_features_filename) and (should_plot_simulation_mean_concentrations or should_plot_variable_mean_concentrations or should_plot_pair_mean_concentrations):
                raise NotImplementedError
        
        
        if not simulations_features_all_included:
            raise RuntimeError('Computed features are not available for some simulations. If running multiple jobs, run this script again when all jobs have completed. If you still get this error, check for and delete .tmp files that may have been left by another run that failed.')
        
        simulations_features = xr.concat(simulations_features, dim=pd.Index(geometry_names, name='Geometry name'))
        
        simulations_iterations_features = xr.concat(simulations_iterations_features, dim=pd.Index(geometry_names, name='Geometry name'))
        
        if should_plot_simulations and len(simulations_iterations_features) > 0:
            molecules_counts = simulations_iterations_features.sel({'Iteration feature': 'molecule_count'}).reset_coords(drop=True)
            molecules_counts.name = 'molecule_count'
            plot_and_save_xarray_combinations(molecules_counts, 'Iteration', ['Geometry name'], '{} plot'.format(data_set_base_filename), (8, 6), should_overwrite=debug_always_overwrite_plots)
        
        
        geometric_coordinates = load_geometric_coordinates(info['mcell_simulations_directory_pattern'])
        geometric_coordinates.coords['Geometry name'] = geometric_coordinates.coords['Geometry name'].pipe(np.vectorize(mcell_geometry_name_abbreviation_func))
        geometric_coordinates = geometric_coordinates.loc[geometry_names, :]
        
        mean_concentrations = []
        
        result = {}
        # result['simulations_base_dir'] = simulations_base_dir
        # result['simulation_data_name_pattern'] = simulation_data_name_pattern
        result['data_set_id'] = data_set_id
        result['data_set_base_filename'] = data_set_base_filename
        result['geometric_coordinates'] = geometric_coordinates
        result['simulations_features'] = simulations_features
        result['selected_variables_names'] = simulations_features.coords['Variable'].values
        result['mean_concentrations'] = mean_concentrations
        result['geometry_names'] = geometry_names
        
        return result
    
    
    def vcell_exported_dir_simulation_id_func(exported_dir):
        return int(re.fullmatch(r'SimID_([0-9]+)_exported.ucd', os.path.split(exported_dir)[1]).group(1))
    
    def vcell_exported_file_token_func(exported_file):
        return re.fullmatch(r'SimID_([0-9]+)_0__vol_([0-9]+).ucd', os.path.split(exported_file)[1]).groups()
    
    def vcell_exported_file_simulation_id_func(exported_file):
        return int(vcell_exported_file_token_func(exported_file)[0])
    
    def vcell_exported_file_iteration_func(exported_file):
        return int(vcell_exported_file_token_func(exported_file)[1])
    
    def vcell_analysis(info, data_set_id):
        '''
        Preprocess simulation data for plotting and regression analysis
        '''
        simulation_data_dir_list = glob.glob(os.path.join(info['vcml_simulations_directory_pattern'], 'SimID_[0-9]*_exported'))
        simulation_data_dir_list = filename_sorted(simulation_data_dir_list)
        simulation_data_file_list = glob.glob(os.path.join(vcml_simulations_directory_pattern, 'SimID_[0-9]*.ucd'))
        simulation_data_file_list = filename_sorted(simulation_data_file_list)
        simulation_jobs_status_filename = os.path.join(vcml_simulations_directory_pattern, 'simulation_jobs_status.txt')
        
        with open(simulation_jobs_status_filename, 'r') as f:
            simulation_jobs_status_filename_contents = f.read()
        raise NotImplementedError('Old dataset-specific path used below')
        simulation_id_to_geometry_name_contents_matches = list(re.finditer(r'^\(BM\) img_.*?, \(App\) (img_SarmaGhosh2012ForCO_0.50_[0-9]{10}_4000sec).*?[Ss]patial.*?\t.+?\t([0-9]+) \/ [0-9]+ \/ [0-9]+\t.+?$', simulation_jobs_status_filename_contents, re.MULTILINE))
        simulation_id_to_geometry_name = {int(x.group(2)): x.group(1) for x in simulation_id_to_geometry_name_contents_matches}
        geometry_name_to_simulation_id = {y: x for x, y in simulation_id_to_geometry_name.items()}
        
        data_set_base_filename = os.path.join(base_filename, data_set_id)
        if debug_max_simulations is not None:
            data_set_base_filename += '_maxsim{0:06d}'.format(debug_max_simulations)
        if debug_max_times is not None:
            data_set_base_filename += '_maxt{0:06d}'.format(debug_max_times)
        if debug_max_faces is not None:
            data_set_base_filename += '_maxf{0:06d}'.format(debug_max_faces)
        
        print('len(simulation_data_dir_list) "{}"'.format(len(simulation_data_dir_list)))
        if debug_max_simulations is not None:
            simulation_data_dir_list = simulation_data_dir_list[:debug_max_simulations]
            print('len(simulation_data_dir_list) "{}"'.format(len(simulation_data_dir_list)))
        
        
        # Get file lists from VCell simulation data exported in UCD format, either in subdirectories or extracted into the same directory
        geometry_names = []
        simulation_data_file_lists = {}
        simulation_base_filenames = {}
        
        for simulation_data_dir in simulation_data_dir_list:
            simulation_data_dir_name = os.path.split(simulation_data_dir)[-1]
            simulation_base_filename = '{} {}'.format(data_set_base_filename, simulation_data_dir_name)
            simulation_id = vcell_exported_dir_simulation_id_func(simulation_data_dir_name)
            
            simulation_features_filename = get_xarray_filename(simulation_base_filename, 'features')
            
            file_path_list = glob.glob(os.path.join(simulation_data_dir, 'SimID_[0-9]*.ucd'))
            file_path_list = filename_sorted(file_path_list)
            simulation_data_file_lists[simulation_id] = file_path_list
            simulation_base_filenames[simulation_id] = simulation_base_filename
        
        print('len(simulation_data_file_lists)', len(simulation_data_file_lists))
        
        for simulation_data_file in simulation_data_file_list:
            simulation_id = vcell_exported_file_simulation_id_func(simulation_data_file)
            simulation_base_filename = '{} {}'.format(data_set_base_filename, simulation_id)
            
            if simulation_id not in simulation_data_file_lists:
                simulation_data_file_lists[simulation_id] = []
            simulation_data_file_lists[simulation_id].append(simulation_data_file)
            simulation_base_filenames[simulation_id] = simulation_base_filename
        
        simulation_ids = sorted(simulation_data_file_lists.keys())
        print('len(simulation_data_file_lists)', len(simulation_data_file_lists))
        
        
        simulations_faces = []
        simulations_vertices = []
        simulations_concentrations_masked = []
        simulations_faces_mask = []
        simulations_concentrations_masked_spatial_means = []
        
        simulations_features = []
        
        for simulation_id in simulation_ids:
            simulation_data_file_list = simulation_data_file_lists[simulation_id]
            simulation_base_filename = simulation_base_filenames[simulation_id]
            geometry_full_name = simulation_id_to_geometry_name[simulation_id]
            geometry_name = mcell_geometry_name_abbreviation_func(geometry_full_name)
            geometry_names.append(geometry_name)
            
            simulation_concentrations_masked_spatial_means_filename = get_xarray_filename(simulation_base_filename, 'concentrations_masked_spatial_means')
            simulation_concentrations_masked_spatial_means_filename = get_xarray_filename(simulation_base_filename, 'concentrations_masked_spatial_means')
            
            simulation_features_filename = get_xarray_filename(simulation_base_filename, 'features')
            
            # Load data
            
            if debug_always_overwrite_data or debug_always_load_features or not os.path.exists(simulation_features_filename):
                if not debug_always_overwrite_data and os.path.exists(simulation_concentrations_masked_spatial_means_filename):
                    simulation_faces = load_xarray(simulation_base_filename, ['faces'])['faces']
                    simulation_vertices = load_xarray(simulation_base_filename, ['vertices'])['vertices']
                    simulation_concentrations_masked = load_xarray(simulation_base_filename, ['concentrations_masked'])['concentrations_masked']
                    simulation_faces_mask = load_xarray(simulation_base_filename, ['faces_mask'])['faces_mask']
                    simulation_concentrations_masked_spatial_means = load_xarray(simulation_base_filename, ['concentrations_masked_spatial_means'])['concentrations_masked_spatial_means']
                    
                else:
                    
                    file_path_list = simulation_data_file_list
                    
                    if debug_max_times is not None:
                        file_path_list = file_path_list[:debug_max_times]
                    
                    is_first_file = True
                    simulation_data_UCD_info = {}
                    variable_names = None
                    for j, file_path in enumerate(file_path_list):
                        include_mesh = j == 0
                        include_mesh = j < 2 # Debug
                        timepoint_UCD_info = parse_UCD_file(file_path, given_include_mesh=include_mesh, given_return_format='xarray')
                        filename = timepoint_UCD_info['filename']
                        variable_names = timepoint_UCD_info['variable_names']
                        simulation_data_UCD_info[filename] = timepoint_UCD_info
                    
                    compute_UCD_preprocessing_result = compute_UCD_preprocessing(simulation_data_UCD_info, given_save_mesh=True, given_ignore_constant_faces=True)
                    
                    simulation_faces = compute_UCD_preprocessing_result['faces']
                    simulation_vertices = compute_UCD_preprocessing_result['vertices']
                    simulation_concentrations_masked = compute_UCD_preprocessing_result['all_face_data']
                    simulation_faces_mask = compute_UCD_preprocessing_result['faces_mask']
                    simulation_concentrations_masked_spatial_means = simulation_concentrations_masked.mean(dim='Face index')
                    
                    save_xarray(simulation_base_filename, faces=simulation_faces)
                    save_xarray(simulation_base_filename, vertices=simulation_vertices)
                    save_xarray(simulation_base_filename, concentrations_masked=simulation_concentrations_masked)
                    save_xarray(simulation_base_filename, faces_mask=simulation_faces_mask)
                    save_xarray(simulation_base_filename, concentrations_masked_spatial_means=simulation_concentrations_masked_spatial_means)
            
            
            # Compute features
            
            if os.path.exists(simulation_features_filename):
                simulation_features = load_xarray(simulation_base_filename, ['features'])['features']
                
            else:
                
                file_path_list = simulation_data_file_list
                
                if debug_max_times is not None:
                    file_path_list = file_path_list[:debug_max_times]
                
                if not debug_skip_computing_features:
                    simulation_features = compute_UCD_simulation_features_20200203(simulation_concentrations_masked)
                
                save_xarray(simulation_base_filename, features=simulation_features)
            
            simulations_features.append(simulation_features)
            
            # Optionally load other data
            
            if os.path.exists(simulation_features_filename) and (should_plot_simulation_mean_concentrations or should_plot_variable_mean_concentrations or should_plot_pair_mean_concentrations):
                simulation_concentrations_masked_spatial_means = load_xarray(simulation_base_filename, ['concentrations_masked_spatial_means'])['concentrations_masked_spatial_means']
                simulations_concentrations_masked_spatial_means.append(simulation_concentrations_masked_spatial_means)
        
        
        
        if len(simulations_concentrations_masked_spatial_means) > 0:
            simulations_concentrations_masked_spatial_means = xr.concat(simulations_concentrations_masked_spatial_means, dim=pd.Index(geometry_names, name='Geometry name'))
        mean_concentrations = simulations_concentrations_masked_spatial_means
        simulations_info = dict(
            simulations_faces=simulations_faces,
            simulations_vertices=simulations_vertices,
            simulations_concentrations_masked=simulations_concentrations_masked,
            simulations_faces_mask=simulations_faces_mask,
            simulations_concentrations_masked_spatial_means=simulations_concentrations_masked_spatial_means,
            simulations_features=simulations_features,
            )
        
        simulations_features = xr.concat(simulations_features, dim=pd.Index(geometry_names, name='Geometry name'))
        
        geometric_coordinates = load_geometric_coordinates(info['vcml_simulations_directory_pattern'])
        geometric_coordinates.coords['Geometry name'] = geometric_coordinates.coords['Geometry name'].pipe(np.vectorize(mcell_geometry_name_abbreviation_func))
        geometric_coordinates = geometric_coordinates.loc[geometry_names, :]
            
        result = {}
        # result['simulations_base_dir'] = simulations_base_dir
        # result['simulation_data_name_pattern'] = simulation_data_name_pattern
        result['data_set_id'] = data_set_id
        result['data_set_base_filename'] = data_set_base_filename
        result['geometric_coordinates'] = geometric_coordinates
        result['simulations_features'] = simulations_features
        # result['selected_variables_names'] = variable_names
        # result['selected_variables_names'] = ['mind_m', 'mind_adp', 'mine', 'minde_m', 'mind_atp']
        result['selected_variables_names'] = simulation_concentrations_masked.coords['Variable'].values
        result['mean_concentrations'] = mean_concentrations
        result['geometry_names'] = geometry_names
        
        return result
    
    
    
    
    
    
    mcell_analysis_load_func = lambda: mcell_analysis(info=info, data_set_id=mcell_data_set_id)
    
    vcell_analysis_load_func = lambda: vcell_analysis(info=info, data_set_id=vcell_data_set_id)
    
    
    # Feature computation
    
    
    analysis_load_funcs = {}
    if enable_mcell_analysis:
        analysis_load_funcs[mcell_data_set_id] = mcell_analysis_load_func
    
    if enable_vcell_analysis:
        raise ValueError('Paths not set')
        analysis_load_funcs[vcell_data_set_id] = vcell_analysis_load_func
    
    
    # Plotting and statistical analysis
    
    figure_base_height = 3
    figure_aspect = 1.6
    figure_base_width = figure_base_height * figure_aspect
    
    default_font_size = 16
    default_marker_size = 12
    default_line_width = 4
    matplotlib.rc('font', size=default_font_size)
    matplotlib.rc('lines', markersize=default_marker_size)
    matplotlib.rc('lines', linewidth=default_line_width)
    
    for data_set_id in sorted(analysis_load_funcs.keys()):
        analysis_load_func = analysis_load_funcs[data_set_id]
        data_set_base_filename = os.path.join(base_filename, data_set_id)
        
        load_result = analysis_load_func()
        data_set_id = load_result['data_set_id']
        data_set_base_filename = load_result['data_set_base_filename']
        geometric_coordinates = load_result['geometric_coordinates']
        simulations_features = load_result['simulations_features']
        selected_variables_names = load_result['selected_variables_names']
        mean_concentrations = load_result['mean_concentrations']
        geometry_names = load_result['geometry_names']
        
        geometric_feature_names = geometric_coordinates.coords[geometric_feature_prefix].values
        
        analysis_filename = f'{data_set_base_filename} analysis'
        can_start, final_exists, temp_name = chunk_start(analysis_filename)
        # if debug_always_overwrite_analysis and final_exists:
        if debug_always_overwrite_analysis:
            can_start = True
        
        if not can_start:
            continue
        
        
        figure_height = figure_base_height * len(geometric_coordinates.coords[geometric_feature_prefix])
        figure_width = figure_height * figure_aspect
        
        
        raw_features = simulations_features
        
        variables_names = raw_features.coords['Variable']
        features_names = raw_features.coords['Feature']
        
        
        control_method_and_feature_names_combinations = []
        if control_method != 'none':
            control_method_and_feature_names_combinations.append((control_method, 'all', geometric_feature_names))
        control_method_and_feature_names_combinations.append(('none', 'all', geometric_feature_names))
        control_method_and_feature_names_combinations.append(('none', 'shapespace', space_space_dimension_names_func(8)))
        control_method_and_feature_names_combinations.append(('none', 'confounders', independent_features_to_control_for))
        
        control_method_backup = str(control_method)
        geometric_coordinates_backup = geometric_coordinates.copy()
        
        for control_method, geometric_feature_names2_name, geometric_feature_names2 in control_method_and_feature_names_combinations:
            
            features_filename_main = data_set_base_filename + ' features'
            features_filename = features_filename_main + '.npz'
            stats_base_filename = features_filename_main + ' stats'
            stats_base_filename += f'_ctrl-{control_method}'
            stats_base_filename += f'_geomfeatset-{geometric_feature_names2_name}'
            stats_reg_base_filename = stats_base_filename + '_reg'
            
            geometric_feature_names3 = list(filter(lambda x: x in geometric_coordinates_backup.coords[geometric_feature_prefix].values, geometric_feature_names2))
            
            
            if len(geometric_feature_names3) > 0:
                geometric_coordinates = geometric_coordinates_backup.loc[{geometric_feature_prefix: geometric_feature_names3}].copy()
                
                features_names = features_names.copy()
                variables_names = variables_names.copy()
                
                should_select_features = False
                # should_select_features = True
                if should_select_features:
                    selected_features_names = features_names.reindex(['mean mean_crossing', 'mean slope'])[0]
                    raw_features = raw_features.loc[{'Feature': selected_features_names}]
                
                # should_select_variables = False
                should_select_variables = True
                if should_select_variables:
                    selected_variables_names = variables_names.loc[selected_variables_names]
                    raw_features = raw_features.loc[{'Variable': selected_variables_names}]
                    if len(mean_concentrations) > 0:
                        mean_concentrations = mean_concentrations.loc[{'Variable': selected_variables_names}]
                
                geometric_coordinates = geometric_coordinates.loc[{'Geometry name': geometry_names}]
                geometric_coordinates_backup = geometric_coordinates.copy()
                geometry_features_notnull = geometric_coordinates.notnull().all(dim='Geometry name')
                print('geometry_features_notnull.sum()', geometry_features_notnull.sum())
                geometric_coordinates = geometric_coordinates.loc[{geometric_feature_prefix: geometry_features_notnull}]
                # debug_raise(analysis_filename, can_start)
                
                
                # Control geometric features, which will be independent variables in the regression analysis
                
                # Only control for a subset of features
                coords = geometric_coordinates.coords[geometric_feature_prefix].values
                # print('coords'); print(coords) # Debug
                # coords_mask = [x for x in coords if any([fnmatch.fnmatchcase(x, y) for y in independent_features_to_control_for])]
                coords_mask = np.array([any([fnmatch.fnmatchcase(x, y) for y in independent_features_to_control_for]) for x in coords])
                coords_with_measurable_mean_mask = np.array([any([fnmatch.fnmatchcase(x, y) for y in independent_features_with_measurable_mean]) for x in coords])
                coords_orth = coords[coords_mask]
                coords_orth_with_measurable_mean = coords[coords_mask & coords_with_measurable_mean_mask]
                geometric_coordinates_to_control_for = []
                for exponent in range(1, independent_features_control_max_exponent + 1):
                    temp = geometric_coordinates.loc[{geometric_feature_prefix: coords_orth}].copy()
                    # Remove mean
                    if len(coords_orth_with_measurable_mean) > 0:
                        temp.loc[{geometric_feature_prefix: coords_orth_with_measurable_mean}] -= temp.loc[{geometric_feature_prefix: coords_orth_with_measurable_mean}].mean(geometric_feature_prefix)
                    temp = temp.assign_coords({geometric_feature_prefix: np.array([f'{x} pow{exponent}' for x in temp.coords[geometric_feature_prefix].values])})
                    temp.values **= exponent
                    geometric_coordinates_to_control_for.append(temp)
                geometric_coordinates_to_control_for = xr.concat(geometric_coordinates_to_control_for, dim=geometric_feature_prefix)
                geometric_coordinates_controlled = geometric_coordinates.copy().loc[{geometric_feature_prefix: coords[~coords_mask]}]
                
                if control_method == 'qr':
                    # Use a QR decomposition to orthogonalize confounding variables and project remaining variables into the null space
                    
                    # geometric_coordinates_to_control_for_q, geometric_coordinates_to_control_for_r = np.linalg.qr(geometric_coordinates_to_control_for.values.T)
                    geometric_coordinates_to_control_for_q, geometric_coordinates_to_control_for_r = np.linalg.qr(geometric_coordinates_to_control_for.values)
                    # print('geometric_coordinates.shape', geometric_coordinates.shape) # Debug
                    # print('geometric_coordinates_to_control_for.shape', geometric_coordinates_to_control_for.shape) # Debug
                    # print('geometric_coordinates_to_control_for_q.shape', geometric_coordinates_to_control_for_q.shape) # Debug
                    # print('geometric_coordinates_to_control_for_r.shape', geometric_coordinates_to_control_for_r.shape) # Debug
                    m, n = geometric_coordinates_to_control_for.shape
                    n2, m2 = geometric_coordinates_to_control_for_q.shape
                    # basis = np.concatenate([geometric_coordinates_to_control_for_q, np.zeros((n - m2, m2))], axis=0)
                    basis = geometric_coordinates_to_control_for_q.copy()
                    # print('m, n, m2, n2'); print(m, n, m2, n2) # Debug
                    # print('basis.shape'); print(basis.shape) # Debug
                    # print('(np.abs(basis) > 1e-10) * 1'); print((np.abs(basis) > 1e-10) * 1) # Debug
                    # print('basis'); print(basis)
                    
                    # Check orthogonality
                    # print('Check orthogonality of geometric_coordinates_to_control_for and geometric_coordinates_to_control_for'); print((np.abs((geometric_coordinates_to_control_for.values.T @ geometric_coordinates_to_control_for_q)) > 1e-10) * 1)
                    
                    # Rejection (removal of projection) of features from basis
                    basis_feature_dot_products = basis.T @ geometric_coordinates.values
                    # print('basis_feature_dot_products.shape'); print(basis_feature_dot_products.shape) # Debug
                    geometric_coordinates_basis_projection = basis @ basis_feature_dot_products
                    # print('geometric_coordinates_basis_projection.shape'); print(geometric_coordinates_basis_projection.shape) # Debug
                    geometric_coordinates_controlled -= geometric_coordinates_basis_projection
                    
                    # Rename features to indicate modification
                    coords2 = coords.copy()
                    coords2[coords_mask] = [f'Or{min(x, geometric_coordinates_to_control_for_q.shape[1])}Pow{independent_features_control_max_exponent} {y}' for x, y in enumerate(coords_orth)]
                    coords2[~coords_mask] = [f'Or{geometric_coordinates_to_control_for_q.shape[1]} {x}' for x in coords[~coords_mask]]
                    geometric_coordinates_controlled.coords[geometric_feature_prefix] = coords2
                    # print('geometric_coordinates_controlled'); print(geometric_coordinates_controlled) # Debug
                    # print('(np.abs(geometric_coordinates_controlled) > 1e-10) * 1'); print((np.abs(geometric_coordinates_controlled) > 1e-10) * 1) # Debug
                    # Remove confounding variables or comment to test for removed relationship (increases p-values some)
                    geometric_coordinates_controlled = geometric_coordinates_controlled.loc[{geometric_feature_prefix: coords2[~coords_mask]}]
                    # debug_raise(analysis_filename, can_start)
                    
                elif control_method == 'ols':
                    # Compute linear regression models of uncontrolled vs. controlled variables and subtract from uncontrolled leaving only residuals. Unlike the `'qr'` method, this can account for constant offsets.
                    
                    # debug_raise(analysis_filename, can_start)
                    
                    # Linear regression
                    geometric_coordinates_controlled_backup = geometric_coordinates_controlled.copy() # Debug
                    geometric_coordinates_to_control_for2 = [geometric_coordinates_to_control_for]
                    geometric_coordinates_to_control_for2_constant = xr.DataArray(np.ones((len(geometry_names), 1)), dims=('Geometry name', geometric_feature_prefix), coords={'Geometry name': geometry_names, geometric_feature_prefix: ['const']})
                    geometric_coordinates_to_control_for2.append(geometric_coordinates_to_control_for2_constant)
                    geometric_coordinates_to_control_for2.append(geometric_coordinates_to_control_for2_constant)
                    geometric_coordinates_to_control_for2 = xr.concat(geometric_coordinates_to_control_for2, dim=geometric_feature_prefix)
                    for feature_name in geometric_coordinates_controlled.coords[geometric_feature_prefix].values:
                        control_results = sm.OLS(geometric_coordinates_controlled.loc[{geometric_feature_prefix: feature_name}].to_pandas(), geometric_coordinates_to_control_for2.to_pandas()).fit()
                        control_results_params = control_results.params.copy()
                        # debug_raise(analysis_filename, can_start)
                        # Removal of residuals
                        geometric_coordinates_controlled.loc[{geometric_feature_prefix: feature_name}] = control_results.resid
                        # debug_raise(analysis_filename, can_start)
                    
                    # debug_raise(analysis_filename, can_start)
                    
                    # Rename features to indicate modification
                    coords2 = coords[~coords_mask].copy()
                    coords2[:] = [f'Res{coords_mask.sum()}_{independent_features_control_max_exponent} {x}' for x in coords[~coords_mask]]
                    geometric_coordinates_controlled.coords[geometric_feature_prefix] = coords2
                    # Remove confounding variables or comment to test for removed relationship (increases p-values some)
                    # geometric_coordinates_controlled = geometric_coordinates_controlled.loc[{geometric_feature_prefix: coords2[~coords_mask]}]
                    
                    # debug_raise(analysis_filename, can_start)
                
                geometric_coordinates = geometric_coordinates_controlled
            
            
            
            
            # Regression
            
            if should_compute_regression_statistics:
                regression_models_count = 0
                independent_features = [geometric_coordinates.rename({geometric_feature_prefix: predictor_feature_generic_name})]
                predictor_feature_names_without_const = independent_features[0].coords[predictor_feature_generic_name].values.copy()
                independent_features_constant = xr.DataArray(np.ones((len(geometry_names), 1)), dims=('Geometry name', predictor_feature_generic_name), coords={'Geometry name': geometry_names, predictor_feature_generic_name: ['const']})
                independent_features.append(independent_features_constant)
                independent_features = xr.concat(independent_features, dim=predictor_feature_generic_name)
                predictor_feature_names_with_const = independent_features.coords[predictor_feature_generic_name].values.copy()
                predictor_feature_names = independent_features[0].coords[predictor_feature_generic_name].values.copy()
                
                dependent_features = raw_features
                # dependent_features = dependent_features.loc[{'Feature': dependent_features.notnull().mean(dim=('Geometry name','Variable',)) >= 0.5}]
                dependent_features_coords_ordered = [dependent_features.coords[x].values for x in dependent_features.dims]
                dependent_features_dims_ordered_without_obs = list(filter(lambda x: x != 'Geometry name', dependent_features.dims))
                dependent_features_coords_ordered_without_obs = [dependent_features.coords[x].values for x in dependent_features_dims_ordered_without_obs]
                dependent_features_shape_ordered_without_obs = [len(x) for x in dependent_features_coords_ordered_without_obs]
                dependent_features_dims_ordered_without_obs_feat = list(filter(lambda x: x not in ['Geometry name', predictor_feature_generic_name], dependent_features.dims))
                dependent_features_coords_ordered_without_obs_feat = [dependent_features.coords[x].values for x in dependent_features_dims_ordered_without_obs_feat]
                dependent_features_shape_ordered_without_obs_feat = [len(x) for x in dependent_features_coords_ordered_without_obs_feat]
    
            
                hypothesis_test_count = 0
                p_value_threshold = 0.01
            
            
                # Regression of single features onto single shape space dimensions
            
                print()
                print('Regression of single features onto single shape space dimensions')
            
                single_feature_coefficients_names = ['const', predictor_feature_generic_name]
                single_feature_coefficients = xr.DataArray(np.nan + np.zeros((len(predictor_feature_names), *dependent_features_shape_ordered_without_obs, len(single_feature_coefficients_names))), dims=(predictor_feature_generic_name, *dependent_features_dims_ordered_without_obs, 'Regression input name',), coords=[predictor_feature_names, *dependent_features_coords_ordered_without_obs, single_feature_coefficients_names])
            
                single_feature_regression_statistics_names = []
                single_feature_regression_statistics_names.append('F-test p-value')
                single_feature_regression_statistics_names.append('Adjusted R-squared')
                single_feature_regression_statistics_names.append('F-test p-value corrected')
                single_feature_regression_statistics_names.append('F-test null rejected')
            
                single_feature_n_hypothesis_tests = 1
            
                single_feature_regressions_statistics = xr.DataArray(np.nan + np.zeros((len(predictor_feature_names_without_const), *dependent_features_shape_ordered_without_obs, len(single_feature_regression_statistics_names))), dims=(predictor_feature_generic_name, *dependent_features_dims_ordered_without_obs, 'Regression statistic'), coords=[predictor_feature_names_without_const, *dependent_features_coords_ordered_without_obs, single_feature_regression_statistics_names])
            
                hypothesis_test_count += single_feature_regressions_statistics.size / len(single_feature_regression_statistics_names) * single_feature_n_hypothesis_tests
            
                single_feature_using_regression_statistic_func = lambda x: x in single_feature_regressions_statistics.coords['Regression statistic'].values
            
                for predictor_feature, dependent_feature_coord in product(predictor_feature_names_without_const, product(*dependent_features_coords_ordered_without_obs)):
                    dependent_feature_coord2 = dict(zip(dependent_features_dims_ordered_without_obs, dependent_feature_coord))
                    dependent_feature_coord3 = dict(dependent_feature_coord2)
                    dependent_feature_coord3[predictor_feature_generic_name] = predictor_feature
                    
                    column_independent_features = independent_features.loc[{predictor_feature_generic_name: ['const', predictor_feature]}]
                    for i, name in enumerate(column_independent_features.coords[predictor_feature_generic_name].values):
                        column_independent_features.coords[predictor_feature_generic_name].values[i] = predictor_feature_generic_name if name.startswith(predictor_feature_generic_name) else name
                    column_results = sm.OLS(dependent_features.loc[dependent_feature_coord2].reset_coords(drop=True).to_series(), column_independent_features.reset_coords(drop=True).to_pandas()).fit()
                    column_results_params = column_results.params.copy()
                    column_results_params.index.values[column_results_params.index.values == predictor_feature] = predictor_feature_generic_name
                        
                    dependent_feature_coord4 = dict(dependent_feature_coord3)
                    dependent_feature_coord4['Regression input name'] = column_results_params.index
                    single_feature_coefficients.loc[dependent_feature_coord4] = column_results_params
                    
                    dependent_feature_coord4 = dict(dependent_feature_coord3)
                    dependent_feature_coord4['Regression statistic'] = 'F-test p-value'
                    single_feature_regressions_statistics.loc[dependent_feature_coord4] = column_results.f_pvalue
                    dependent_feature_coord4['Regression statistic'] = 'Adjusted R-squared'
                    single_feature_regressions_statistics.loc[dependent_feature_coord4] = column_results.rsquared_adj
                    
                    regression_models_count += 1
                
                temp1 = single_feature_coefficients.copy().drop('const', dim='Regression input name').reset_coords(drop=True)
                temp2 = independent_features.copy().drop('const', dim=predictor_feature_generic_name)
                single_feature_coefficients_relative_std = (temp1 / (dependent_features.std('Geometry name') / temp2.std('Geometry name'))).std()
                print('single_feature_coefficients_relative_std = {0}'.format(float(single_feature_coefficients_relative_std)))
            
            
                # Regression of single features onto all shape space dimensions
            
                print()
                print('Regression of single features onto all shape space dimensions')
            
                all_features_coefficients = xr.DataArray(np.nan + np.zeros((len(predictor_feature_names), *dependent_features_shape_ordered_without_obs)), dims=(predictor_feature_generic_name, *dependent_features_dims_ordered_without_obs), coords=[predictor_feature_names, *dependent_features_coords_ordered_without_obs])
            
                all_features_regression_statistics_names = []
                all_features_regression_statistics_names.append('F-test p-value')
                all_features_regression_statistics_names.append('Adjusted R-squared')
                all_features_regression_statistics_names.append('F-test p-value corrected')
                all_features_regression_statistics_names.append('F-test null rejected')
            
                all_features_n_hypothesis_tests = 1
            
                all_features_regressions_statistics = xr.DataArray(np.nan + np.zeros((*dependent_features_shape_ordered_without_obs_feat, len(single_feature_regression_statistics_names))), dims=(*dependent_features_dims_ordered_without_obs_feat, 'Regression statistic'), coords=[*dependent_features_coords_ordered_without_obs_feat, single_feature_regression_statistics_names])
            
                hypothesis_test_count += all_features_regressions_statistics.size / len(all_features_regression_statistics_names) * all_features_n_hypothesis_tests
            
                all_features_using_regression_statistic_func = lambda x: x in all_features_regressions_statistics.coords['Regression statistic'].values
                
                for dependent_feature_coord in product(*dependent_features_coords_ordered_without_obs_feat):
                    dependent_feature_coord2 = dict(zip(dependent_features_dims_ordered_without_obs_feat, dependent_feature_coord))
                    
                    column_results = sm.OLS(dependent_features.loc[dependent_feature_coord2].reset_coords(drop=True).to_series(), independent_features.to_pandas()).fit()
                    column_results_params = column_results.params.copy()
                        
                    all_features_coefficients.loc[dependent_feature_coord2] = column_results_params
                    
                    dependent_feature_coord3 = dict(dependent_feature_coord2)
                    dependent_feature_coord3['Regression statistic'] = 'F-test p-value'
                    all_features_regressions_statistics.loc[dependent_feature_coord3] = column_results.f_pvalue
                    dependent_feature_coord3['Regression statistic'] = 'Adjusted R-squared'
                    all_features_regressions_statistics.loc[dependent_feature_coord3] = column_results.rsquared_adj
                    
                    regression_models_count += 1
                
                temp1 = all_features_coefficients.copy().drop('const', dim=predictor_feature_generic_name).reset_coords(drop=True)
                temp2 = independent_features.copy().drop('const', dim=predictor_feature_generic_name)
                all_features_coefficients_relative_std = (temp1 / (dependent_features.std('Geometry name') / temp2.std('Geometry name'))).std()
                print('all_features_coefficients_relative_std = {0}'.format(float(all_features_coefficients_relative_std)))
            
                # Print statistics with comprehensive Bonferroni correction
            
                regressions_statistics_base_filename = stats_reg_base_filename
            
                single_feature_f_pvalues = single_feature_regressions_statistics.loc[{'Regression statistic': 'F-test p-value'}]
                single_feature_f_pvalues_multipletests_results = single_feature_f_pvalues.copy()
                single_feature_f_pvalues_corrected_reject = single_feature_f_pvalues.copy()
                single_feature_f_pvalues_corrected = single_feature_f_pvalues.copy()
                single_feature_f_pvalues_corrected *= hypothesis_test_count
                single_feature_f_pvalues_corrected = single_feature_f_pvalues_corrected.clip(max=1)
                single_feature_f_pvalues_corrected_reject = single_feature_f_pvalues_corrected <= p_value_threshold
                single_feature_regressions_statistics.loc[{'Regression statistic': 'F-test p-value corrected'}] = single_feature_f_pvalues_corrected
                single_feature_regressions_statistics.loc[{'Regression statistic': 'F-test null rejected'}] = single_feature_f_pvalues_corrected_reject
            
                all_features_f_pvalues = all_features_regressions_statistics.loc[{'Regression statistic': 'F-test p-value'}]
                all_features_f_pvalues_multipletests_results = all_features_f_pvalues.copy()
                all_features_f_pvalues_corrected_reject = all_features_f_pvalues.copy()
                all_features_f_pvalues_corrected = all_features_f_pvalues.copy()
                all_features_f_pvalues_corrected *= hypothesis_test_count
                all_features_f_pvalues_corrected = all_features_f_pvalues_corrected.clip(max=1)
                all_features_f_pvalues_corrected_reject = all_features_f_pvalues_corrected <= p_value_threshold
                all_features_regressions_statistics.loc[{'Regression statistic': 'F-test p-value corrected'}] = all_features_f_pvalues_corrected
                all_features_regressions_statistics.loc[{'Regression statistic': 'F-test null rejected'}] = all_features_f_pvalues_corrected_reject
                print('regression_models_count = {0}'.format(regression_models_count))
                print('hypothesis_test_count = {0}'.format(hypothesis_test_count))
                
                sort_statistic = 'F-test p-value'
                # sort_statistic = 'F-test p-value corrected'
                
                print()
                print('Saving single_feature_regressions_statistics')
                single_feature_regressions_statistics_stacked = single_feature_regressions_statistics.stack({'Inputs stacked': single_feature_regressions_statistics.dims[:-1]})
                single_feature_regressions_statistics_stacked = single_feature_regressions_statistics_stacked.transpose('Inputs stacked', 'Regression statistic')
                inputs_stacked_multiindex = single_feature_regressions_statistics_stacked.indexes['Inputs stacked']
                inputs_stacked_string_index = pd.Index([', '.join(x) for x in inputs_stacked_multiindex])
                single_feature_regressions_statistics_stacked = single_feature_regressions_statistics_stacked.assign_coords(**{'Inputs stacked': inputs_stacked_string_index})
                single_feature_regressions_statistics_stacked_df = single_feature_regressions_statistics_stacked.to_dataframe(name='')
                single_feature_regressions_statistics_stacked_df.to_csv('{0} single_feature.csv'.format(regressions_statistics_base_filename))
                single_feature_regressions_statistics_stacked_sorted = single_feature_regressions_statistics_stacked.sortby(single_feature_regressions_statistics_stacked.loc[{'Regression statistic': sort_statistic}])
                single_feature_regressions_statistics_df = single_feature_regressions_statistics_stacked.to_dataframe('').unstack()
                single_feature_regressions_statistics_df.columns = single_feature_regressions_statistics_df.columns.levels[1]
                single_feature_regressions_statistics_df_sorted = single_feature_regressions_statistics_df.sort_values(by=[sort_statistic])
                single_feature_regressions_statistics_df_sorted.to_csv('{0} single_feature sorted.csv'.format(regressions_statistics_base_filename))
            
                print()
                print('Saving all_features_regressions_statistics')
                all_features_regressions_statistics_stacked = all_features_regressions_statistics.stack({'Inputs stacked': all_features_regressions_statistics.dims[:-1]})
                all_features_regressions_statistics_stacked = all_features_regressions_statistics_stacked.transpose('Inputs stacked', 'Regression statistic')
                inputs_stacked_multiindex = all_features_regressions_statistics_stacked.indexes['Inputs stacked']
                inputs_stacked_string_index = pd.Index([', '.join(x) for x in inputs_stacked_multiindex])
                all_features_regressions_statistics_stacked = all_features_regressions_statistics_stacked.assign_coords(**{'Inputs stacked': inputs_stacked_string_index})
                all_features_regressions_statistics_stacked_df = all_features_regressions_statistics_stacked.to_dataframe(name='')
                all_features_regressions_statistics_stacked_df.to_csv('{0} all_features.csv'.format(regressions_statistics_base_filename))
                all_features_regressions_statistics_stacked_sorted = all_features_regressions_statistics_stacked.sortby(all_features_regressions_statistics_stacked.loc[{'Regression statistic': sort_statistic}])
                # unstack destroys the sorting
                all_features_regressions_statistics_df = all_features_regressions_statistics_stacked.to_dataframe('').unstack()
                all_features_regressions_statistics_df.columns = all_features_regressions_statistics_df.columns.levels[1]
                all_features_regressions_statistics_df_sorted = all_features_regressions_statistics_df.sort_values(by=[sort_statistic])
                all_features_regressions_statistics_df_sorted.to_csv('{0} all_features sorted.csv'.format(regressions_statistics_base_filename))
            
            
                print()
                print('Saving to netcdf')
                save_xarray(stats_base_filename,
                    single_feature_coefficients=single_feature_coefficients,
                    single_feature_regressions_statistics=single_feature_regressions_statistics,
                    single_feature_regressions_statistics_stacked_sorted=single_feature_regressions_statistics_stacked_sorted,
                    # single_feature_relative_coefficients=single_feature_relative_coefficients,
                    single_feature_coefficients_relative_std=single_feature_coefficients_relative_std,
                    all_features_coefficients=all_features_coefficients,
                    all_features_regressions_statistics=all_features_regressions_statistics,
                    # all_features_relative_coefficients=all_features_relative_coefficients,
                    all_features_coefficients_relative_std=all_features_coefficients_relative_std,
                    all_features_regressions_statistics_stacked_sorted=all_features_regressions_statistics_stacked_sorted,
                    )
            
            
            
            
            # Plot
            
            
            # Scatter plots with coloring
            
            Y_trim = 2.5e-2
            
            print('Y_trim {0}'.format(Y_trim))
            
            
            geometric_feature_prefix2 = 'GF dim.'
            pair_grid_max_chars_per_line = 7
            def name_display_transform_func(given_name, given_max_chars_per_line=None):
                result = str(given_name)
                if result.startswith(geometric_feature_prefix):
                    # result = os.linesep.join(['Dim', result.split(geometric_feature_prefix)[1]])
                    # result = '{0} {1}'.format(geometric_feature_prefix2, result.split(geometric_feature_prefix)[1].trim())
                    result = result.replace(geometric_feature_prefix, geometric_feature_prefix2)
                if given_max_chars_per_line is not None:
                    result = [result]
                    while len(result[-1]) > given_max_chars_per_line:
                        result.append(result[-1][given_max_chars_per_line:])
                        result[-2] = result[-2][:given_max_chars_per_line]
                    result = '\n'.join(result)
                return result
            
            
            
            plot_extensions = ['png']
            # plot_extensions = ['png', 'svg']
            
            # Plot pairs of independent and dependent features
            
            if should_plot_pairs:
                limits_padding = matplotlib.rcParamsDefault['lines.markersize'] / 72
                independent_features_names_to_plot = independent_features.coords[predictor_feature_generic_name].values
                independent_features_names_to_plot = independent_features_names_to_plot[independent_features_names_to_plot != 'const']
                for variable_name, independent_feature_name, dependent_feature_name in product(dependent_features.coords['Variable'].values, independent_features_names_to_plot, dependent_features.coords['Feature'].values):
                    X = independent_features.loc[{predictor_feature_generic_name: independent_feature_name}]
                    Y = dependent_features.loc[{'Variable': variable_name, 'Feature': dependent_feature_name}]
                    
                    figure_title = 'Variable {0}'.format(variable_name)
                    
                    pair_plot_filenames = ['{0} {1} {2} {3}.{4}'.format(features_filename_main, variable_name, independent_feature_name, dependent_feature_name, extension) for extension in plot_extensions]
                    should_plot_formats = [not os.path.exists(x) for x in pair_plot_filenames]
                    
                    if any(should_plot_formats) or debug_always_overwrite_plots:
                        
                        height = 6
                        aspect = 1.6
                        fig_width = height * aspect + default_font_size / 72 * 7
                        fig_height = height + default_font_size / 72 * 5
                        fig = plt.figure(figsize=(fig_width, fig_height))
                        fig_axes = plt.scatter(X, Y)
                        plt.xlabel(independent_feature_name)
                        plt.ylabel(dependent_feature_name)
                        axes_width = fig_axes.axes.get_position().width
                        axes_height = fig_axes.axes.get_position().height
                        xlim_padding = (X.max() - X.min()) * (limits_padding / axes_width)
                        ylim_padding = (Y.max() - Y.min()) * (limits_padding / axes_height)
                        plt.xlim(X.min() - xlim_padding, X.max() + xlim_padding)
                        plt.ylim(Y.min() - ylim_padding, Y.max() + ylim_padding)
                        plt.tight_layout()
                
                        plt.title(figure_title)
                        
                        for pair_plot_filename in pair_plot_filenames:
                            plt.savefig(pair_plot_filename)
                        
                        del fig, fig_axes
                        plt.close('all')
            
            
            # Plot pair grids of shape space coordinates colored by features
            
            palette = sns.mpl_palette('viridis', n_colors=1000)
            cmap = matplotlib.colors.ListedColormap(palette)
            
            if should_plot_pair_grids:
                X = independent_features.drop('const', dim=predictor_feature_generic_name)
                predictor_features_names = X.coords[predictor_feature_generic_name].values
                for variable_name, dependent_feature_name in product(dependent_features.coords['Variable'].values, dependent_features.coords['Feature'].values):
                    Y = dependent_features.loc[{'Variable': variable_name, 'Feature': dependent_feature_name}]
            
                    Y_limits = (Y.quantile(Y_trim), Y.quantile(1 - Y_trim))
                    
                    # Attempt to prevent matplotlib savefig error "minvalue must be less than or equal to maxvalue"
                    for i in range(1, 6):
                        if Y_limits[0] < Y_limits[1]:
                            break
                        Y_trim2 = Y_trim * 2**(-i) if i < 5 else 0
                        Y_limits = (Y.quantile(Y_trim2), Y.quantile(1 - Y_trim2))
                    if Y_limits[0] >= Y_limits[1]:
                        Y_limits[0] = np.median() - 0.5
                        Y_limits[1] = np.median() + 0.5
            
                    grid_plot_filenames = ['{0} {1} {2} pairgrid.{3}'.format(features_filename_main, variable_name, dependent_feature_name, extension) for extension in plot_extensions]
                    should_plot_formats = [not os.path.exists(x) for x in grid_plot_filenames]
                    
                    if any(should_plot_formats) or debug_always_overwrite_plots:
                        n_X_columns = len(X.coords[predictor_feature_generic_name])
                        height = 3
                        aspect = 1.6
                        fig = plt.figure(figsize=(height * aspect * (n_X_columns - 1), height * (n_X_columns - 1)))
                        axis_grid = {x: {x: None for x in predictor_features_names} for x in predictor_features_names}
                        for i1, c1 in enumerate(predictor_features_names):
                            for i2, c2 in enumerate(predictor_features_names):
                                if i2 >= i1:
                                    continue
                                fig_subplot = plt.subplot(n_X_columns - 1, n_X_columns - 1, (i1 - 1) * (n_X_columns - 1) + i2 + 1)
                                axis_grid[c2][c1] = fig_subplot.scatter(X.loc[:, c2], X.loc[:, c1], c=Y, vmin=Y_limits[0], vmax=Y_limits[1], cmap=cmap)
                                fig_subplot.set_frame_on(False)
                                fig_subplot.set_xticks([])
                                fig_subplot.set_yticks([])
                                if i1 == n_X_columns - 1:
                                    plt.xlabel(name_display_transform_func(c2, given_max_chars_per_line=pair_grid_max_chars_per_line))
                                if i2 == 0:
                                    plt.ylabel(name_display_transform_func(c1, given_max_chars_per_line=pair_grid_max_chars_per_line))
                
                        # https://stackoverflow.com/questions/44641669/scatterplot-with-point-colors-representing-a-continuous-variable-in-seaborn-face
                        fig.subplots_adjust(right=0.92)
                        colorbar_axes = fig.add_axes([0.94, 0.2, 0.02, 0.6])
                        figure_title = 'Variable {0}, feature {1}'.format(variable_name, dependent_feature_name)
                        fig.suptitle(figure_title)
                        colorbar_points = plt.scatter([], [], c=[], vmin=Y_limits[0], vmax=Y_limits[1], cmap=cmap)
                        fig.colorbar(colorbar_points, cax=colorbar_axes)
                        for grid_plot_filenames in grid_plot_filenames:
                            plt.savefig(grid_plot_filename)
                        del fig, fig_subplot, axis_grid, colorbar_axes, colorbar_points
                        plt.close('all')
    
        chunk_finish(analysis_filename)
    
    
    
    # debug_raise() # Debug


if __name__ == '__main__':
    
    # Parse options
    substitution_key_set = bsf.cellorganizer_substitution_key_set
    path_with_substitution = bsf.cellorganizer_path_with_substitution
    paths_with_substitution = bsf.cellorganizer_path_with_substitution
    
    parser = argparse.ArgumentParser(description='Analyze the results of MCell or Virtual Cell simulations run by `generate_and_simulate.py`.')
    
    # parser.add_argument('reaction_network_file', type=existing_glob, help='Either a path of a single VCML file or a glob pattern for a single collection of MCell MDL files. Same value given to `generate_and_simulate.py`.')
    parser.add_argument('reaction_network_file', type=path_with_substitution, help='Either a path of a single VCML file or a glob pattern for a single collection of MCell MDL files. Same value given to `generate_and_simulate.py`.')
    
    parser.add_argument('output_dir', type=pathlib.Path, help='The path of a directory into which intermediate and final results will be written Same value given to `generate_and_simulate.py`.')
    
    parser.add_argument('--overwrite', action='store_true', help='Overwrite existing intermediate files instead of failing with error (not implemented)')
    
    parser.add_argument('--dry_run', action='store_true', help='Print actions instead of taking them')
    
    parser.add_argument('--cellorganizer', type=existing_path, default='..', help="The path of a CellOrganizer installation")
    
    parser.add_argument('--synthesis', type=str, choices=['framework', 'all'], default='framework', help="`'all'` to synthesize framework and vesicles or `'framework'` to ignore the vesicle model")
    
    parser.add_argument('--downsampling', type=positive_float, default=1/2, help='Downsampling factor to trade disk and computational requirements for geometric precision')
    
    parser.add_argument('--vcml_relative_downsampling', type=positive_float, default=1/2, help='Factor for further downsampling for writing VCML files. VCell simulation data from fewer than 100 simulations at full resolution can fill available disk space on the VCell server, but this is less of a concern when running simulations on your own cluster.')
    
    parser.add_argument('--simulation_end_time', type=positive_float, default=4000, help='Time at which to end simulation in seconds')
    
    argv = sys.argv
    argv = remove_prepended_arguments(argv)
    args = parser.parse_args(argv[1:])
    args_vars = vars(args)
    
    bsf.cellorganizer_path_substitutions(parser, args_vars)
    
    
    generate_and_simulate_analysis(**args_vars)
    
