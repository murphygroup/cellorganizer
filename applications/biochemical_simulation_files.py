'''
Python module for analyzing simulation results.

Taraz Buck
2019-10-02
'''

import os
import sys
import glob
import re
from math import *
from itertools import product
import struct
import pprint

import numpy as np
import pandas as pd
import xarray as xr
import scipy as sp
import scipy.signal, scipy.fftpack
# import trimesh
from sklearn.preprocessing import StandardScaler
# from sklearn.decomposition import PCA
from sklearn.linear_model import LinearRegression, TheilSenRegressor
# from sklearn.model_selection import RepeatedKFold, cross_val_score
from sklearn.model_selection import RepeatedKFold
from sklearn.metrics import r2_score
from sklearn.pipeline import Pipeline
from sklearn_pandas import cross_val_score
# import statsmodels.api as sm

import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import seaborn as sns

import inspect
# base_filename = os.path.splitext(inspect.stack()[-1][1])[0]
base_filename = os.path.splitext(__file__)[0]
try:
    os.mkdir(base_filename)
except OSError:
    pass



is_pandas = lambda x: isinstance(x, pd.Series) or isinstance(x, pd.DataFrame) or isinstance(x, pd.Panel)
is_xarray = lambda x: isinstance(x, xr.DataArray) or isinstance(x, xr.Dataset)
should_convert_to_numpy_func = lambda x: not is_pandas(x) and not is_xarray(x)

all_but_axis_slice_func = lambda ndim, axis: (slice(None),) * axis + (None,) + (slice(None),) * (ndim - axis - 1)
vector_to_broadcast_slice_func = lambda ndim, axis: (None,) * axis + (slice(None),) + (None,) * (ndim - axis - 1)
"""
all_but_axis_slice_func = lambda ndim, axis: sum([(None,) if (x == axis or (isinstance(x, tuple) and x in axis)) else (slice(None),) for x in range(ndim)], ())
vector_to_broadcast_slice_func = lambda ndim, axis: sum([(slice(None),) if (x == axis or (isinstance(x, tuple) and x in axis)) else (None,) for x in range(ndim)], ())
"""

def mean_func(given_value, given_axis=None):
    if should_convert_to_numpy_func(given_value):
        given_value = np.array(given_value)
    if is_xarray(given_value):
        return_value = given_value.mean(dim=given_axis)
    else:
        return_value = given_value.mean(axis=given_axis)
    return return_value

def std_func(given_value, given_axis=None):
    if should_convert_to_numpy_func(given_value):
        given_value = np.array(given_value)
    if is_xarray(given_value):
        return_value = given_value.std(dim=given_axis)
    else:
        return_value = given_value.std(axis=given_axis)
    return return_value

def contrast_func(given_value, given_quantile, given_axis=None):
    if should_convert_to_numpy_func(given_value):
        given_value = np.array(given_value)
        return_value = np.quantile(given_value, 1 - given_quantile, axis=given_axis) - np.quantile(given_value, given_quantile, axis=given_axis)
    elif is_xarray(given_value):
        return_value = given_value.quantile(1 - given_quantile, dim=given_axis) - given_value.quantile(given_quantile, dim=given_axis)
    else:
        return_value = given_value.quantile(1 - given_quantile, axis=given_axis) - given_value.quantile(given_quantile, axis=given_axis)
    return return_value

def slope_func(given_value, given_edge_order=None, given_axis=None):
    if should_convert_to_numpy_func(given_value):
        given_value = np.array(given_value)
        return_value = np.gradient(given_value, edge_order=2, axis=given_axis)
    elif is_xarray(given_value):
        return_value = given_value.differentiate(coord=given_axis, edge_order=2)
    else:
        return_value = given_value.copy()
        return_value[:] = np.gradient(given_value, edge_order=2, axis=given_axis)
        pass
    return return_value

def mean_crossing_func(given_value):
    given_value_original = given_value
    given_value = np.array(given_value)
    # given_value_mean = given_value.mean()
    given_value_mean = given_value.mean(axis=0)
    return_value = np.zeros(given_value.shape, dtype='bool')
    if return_value.size > 1:
        # return_value[1:] = (given_value[1:] >= given_value_mean) != (given_value[:-1] >= given_value_mean)
        return_value[1:, :] = (given_value[1:, :] >= given_value_mean) != (given_value[:-1, :] >= given_value_mean[None, :])
    if not should_convert_to_numpy_func(given_value_original):
        return_value2 = given_value_original.copy()
        return_value2[:] = return_value
        return_value = return_value2
    return return_value

def mean_above_func(given_value):
    given_value_original = given_value
    given_value = np.array(given_value)
    # given_value_mean = given_value.mean()
    given_value_mean = given_value.mean(axis=0)
    return_value = given_value >= given_value_mean
    if not should_convert_to_numpy_func(given_value_original):
        return_value2 = given_value_original.copy()
        return_value2[:] = return_value
        return_value = return_value2
    return return_value

def corr_func(given_value1, given_value2, given_axis=None):
    if is_xarray(given_value1) and is_xarray(given_value2):
        # Olkin and Pratt 1958, "Unbiased Estimation of Certain Correlation Coefficients"
        """
        given_axis_index1 = (np.array(given_value1.dims) == given_axis).nonzero()[0][0]
        # given_axis_index2 = (np.array(given_value2.dims) == given_axis).nonzero()[0][0]
        given_value1_normalized = (given_value1 - given_value1.mean(dim=given_axis)).values
        given_value2_normalized = (given_value2 - given_value2.mean(dim=given_axis)).values
        given_value1_std = given_value1.std(dim=given_axis).values
        given_value2_std = given_value2.std(dim=given_axis).values
        given_values_mask = (given_value1_std > 0) & (given_value2_std > 0)
        given_value1_broadcast = all_but_axis_slice_func(given_value1.ndim, given_axis_index1)
        given_value1_normalized = given_value1_normalized / given_value1_std[given_value1_broadcast]
        given_value2_normalized = given_value2_normalized / given_value2_std[given_value1_broadcast]
        given_value_correlation = given_value1_normalized * given_value2_normalized
        given_value_correlation_n = np.isfinite(given_value_correlation).sum()
        given_value_correlation = given_value_correlation.sum(axis=given_axis_index1)
        given_value_correlation_mask = np.isfinite(given_value_correlation)
        given_value_correlation2 = given_value_correlation.copy()
        given_value_correlation2[given_value_correlation_mask] *= sp.special.hyp2f1(1/2, 1/2, (int(given_value_correlation_n) - 1) / 2, 1 - given_value_correlation2[given_value_correlation_mask]**2)
        return_value = given_value_correlation
        """
        
        """
        # Olkin and Pratt 1958, "Unbiased Estimation of Certain Correlation Coefficients"
        given_value1_normalized = given_value1 - given_value1.mean(dim=given_axis)
        given_value2_normalized = given_value2 - given_value2.mean(dim=given_axis)
        # given_value1_normalized = given_value1.copy()
        # given_value2_normalized = given_value2.copy()
        given_value1_std = given_value1_normalized.std(dim=given_axis)
        given_value2_std = given_value2_normalized.std(dim=given_axis)
        given_value1_std = given_value1_std.where(given_value1_std > 0, np.nan)
        given_value2_std = given_value2_std.where(given_value2_std > 0, np.nan)
        given_value1_normalized = given_value1_normalized.where(given_value1_std.notnull(), np.nan)
        given_value2_normalized = given_value2_normalized.where(given_value2_std.notnull(), np.nan)
        given_value1_normalized /= given_value1_std
        given_value2_normalized /= given_value2_std
        given_value1_normalized -= given_value1_normalized.mean(dim=given_axis)
        given_value2_normalized -= given_value2_normalized.mean(dim=given_axis)
        given_value_correlation = given_value1_normalized * given_value2_normalized
        given_value_correlation_n = given_value_correlation.notnull().sum()
        given_value_correlation = given_value_correlation.sum(dim=given_axis)
        # given_value_correlation = given_value_correlation * xr.apply_ufunc(lambda x: sp.special.hyp2f1(1/2, 1/2, (given_value_correlation_n - 1) / 2, 1 - x**2), given_value_correlation, input_core_dims=[given_value_correlation.dims])
        given_value_correlation_mask = given_value_correlation.notnull().values
        given_value_correlation_values = given_value_correlation.values.copy()
        given_value_correlation_values[given_value_correlation_mask] *= sp.special.hyp2f1(1/2, 1/2, (int(given_value_correlation_n) - 1) / 2, 1 - given_value_correlation_values[given_value_correlation_mask]**2)
        return_value = given_value_correlation
        """
        
        # Olkin and Pratt 1958, "Unbiased Estimation of Certain Correlation Coefficients"
        given_value1_centered = given_value1 - given_value1.mean(dim=given_axis)
        given_value2_centered = given_value2 - given_value2.mean(dim=given_axis)
        given_values_numerator = (given_value1_centered * given_value2_centered).sum(dim=given_axis)
        given_values_denominator = xr.ufuncs.sqrt(((given_value1_centered**2).sum(dim=given_axis) * (given_value2_centered**2).sum(dim=given_axis)))
        given_values_correlation = given_values_numerator / given_values_denominator
        given_values_correlation = given_values_correlation.where(xr.ufuncs.isfinite(given_values_correlation), np.nan)
        given_values_correlation_n = (xr.ufuncs.isfinite(given_value1_centered) * xr.ufuncs.isfinite(given_value2_centered)).sum(dim=given_axis)
        given_values_correlation_n = given_values_correlation_n.where(xr.ufuncs.isfinite(given_values_correlation), np.nan)
        given_values_correlation2 = given_values_correlation * sp.special.hyp2f1(1/2, 1/2, (given_values_correlation_n - 1) / 2, 1 - given_values_correlation**2)
        given_values_correlation2 = given_values_correlation2.where(xr.ufuncs.isfinite(given_values_correlation2), np.nan)
        return_value = given_values_correlation2
    else:
        raise NotImplementedError
    return return_value


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

cellorganizer_substitution_key_set = {'cellorganizer', 'models', 'data', 'applications'}
def cellorganizer_substitutions_dict_func(cellorganizer_path):
    result = {}
    result['cellorganizer'] = cellorganizer_path
    for subdir in cellorganizer_substitution_key_set - {'cellorganizer'}:
        result[subdir] = os.path.join(result['cellorganizer'], subdir)
    return result

def path_with_substitution(value, substitution_key_set=set()):
    placeholders_pattern = '\{(' + '|'.join(substitution_key_set) + ')\}'
    x_placeholders = re.findall(placeholders_pattern, value)
    if not set(x_placeholders).issubset(substitution_key_set):
        placeholders_doc = '{' + ', '.join(['\'{' + x + '}\'' for x in sorted(substitution_key_set)]) + '}'
        raise argparse.ArgumentTypeError(f'Argument must be an existing path or glob pattern matching at least one existing path with optional placeholders {placeholders_doc} indicating paths relative to the CellOrganizer installation (given {repr(value)})')
    # value = value.format(**substitutions_dict)
    # value = existing_path(value)
    return value

cellorganizer_path_with_substitution = lambda x: path_with_substitution(x, cellorganizer_substitution_key_set)
cellorganizer_paths_with_substitution = lambda x: [path_with_substitution(y, cellorganizer_substitution_key_set) for y in x]

def path_substitution_dict(substitutions_dict, substitution_key_set=set()):
    if not isinstance(substitution_key_set, set):
        substitution_key_set = set(substitution_key_set)
    substitutions_dict2 = {}
    for x in substitution_key_set:
        if x in substitutions_dict:
            substitutions_dict2[x] = substitutions_dict[x]
    return substitutions_dict2

def path_substitution(value, substitutions_dict, substitution_key_set=set()):
    substitutions_dict2 = path_substitution_dict(substitutions_dict, substitution_key_set)
    if isinstance(value, list):
        value = [x.format(**substitutions_dict2) for x in value]
    else:
        value = value.format(**substitutions_dict2)
    # value = existing_path(value)
    return value

def cellorganizer_path_substitutions(parser, args_vars):
    '''
    Performs path substitution on args_vars in place where parser._actions[*].type is in (cellorganizer_path_with_substitution, cellorganizer_paths_with_substitution)
    '''
    substitution_key_set = cellorganizer_substitution_key_set
    substitutions_dict = cellorganizer_substitutions_dict_func(str(args_vars['cellorganizer']))
    
    substitutions_dict2 = path_substitution_dict(substitutions_dict, substitution_key_set)
    
    args_to_substitute = [x.dest for x in parser._actions if x.type in (cellorganizer_path_with_substitution, cellorganizer_paths_with_substitution)]
    # args_vars = {x: y for x, y in args_vars.items()}
    for arg in args_to_substitute:
        value = args_vars[arg]
        if isinstance(value, list):
            value = [x.format(**substitutions_dict2) for x in value]
        else:
            value = value.format(**substitutions_dict2)
        args_vars[arg] = value
    # return args_vars


def filename_sort_key(given_string):
    '''
    Splits on path separators, then whitespace, then finds floats, then splits on underscores and periods.
    '''
    
    given_string = os.path.normpath(given_string)
    fields = list(os.path.split(given_string))
    
    while True:
        temp1, temp2 = os.path.split(fields[0])
        if len(temp1) == 0 or len(temp2) == 0:
            break
        fields = [temp1, temp2] + fields[1:]
    
    fields2 = []
    for field in fields:
        field_fields = field.split()
        fields2.extend(field_fields)
    fields = fields2
    
    fields2 = []
    for field in fields:
        field_fields = []
        while len(field) > 0:
            search_result = float_pattern_program.search(field)
            if search_result is None:
                field_fields.append(field)
                break
            else:
                if search_result.start() > 0:
                    field_fields.append(field[:search_result.start()])
                field_fields.append(search_result.group(0))
                field = field[search_result.end():]
        fields2.extend(field_fields)
    fields = fields2
    
    delimiters = ('_', '.')
    
    for delimiter in delimiters:
        fields2 = []
        for field in fields:
            if float_pattern_program_matches(field) or field in ['.', '..']:
                fields2.append(field)
            else:
                fields3 = field.split(delimiter)
                fields3 = filter(lambda x: len(x) > 0, fields3)
                fields2.extend(fields3)
        fields = fields2
    
    fields2 = []
    for field in fields:
        if field.isdecimal():
            fields2.append(int(field))
        elif float_pattern_program_matches(field):
            fields2.append(float(field))
        else:
            fields2.append(field)
    fields = fields2
    
    # raise NotImplementedError
    
    return fields


def filename_sorted(given_strings):
    return sorted(given_strings, key=filename_sort_key)


def parse_UCD_file(given_file_path, given_include_mesh=False, given_return_format='numpy'):
    '''
    References:
    
    * https://github.com/virtualcell/vcell/blob/master/vcell-core/src/main/java/cbit/vcell/export/AVS_UCD_Exporter.java, writeUCDMembGeomAndData()
    * https://github.com/virtualcell/vcell/blob/master/vcell-core/src/main/java/cbit/vcell/solvers/CartesianMesh.java
    '''
    
    # given_filename = os.path.splitext(given_file_path)
    given_filename = os.path.split(given_file_path)[1]
    
    with open(given_file_path, 'r') as given_file:
        given_file_text_lines = given_file.readlines()
    given_file_text_lines = map(lambda x: x.strip(), given_file_text_lines)
    given_file_text_lines = filter(lambda x: len(x) > 0, given_file_text_lines)
    given_file_text_lines = list(given_file_text_lines)
    
    result = {}
    result['file_path'] = given_file_path
    result['filename'] = given_filename
    
    result['n_vertices'], result['n_faces'], result['n_vertex_data'], result['n_face_data'], result['n_model_data'] = map(int, given_file_text_lines[0].split())
    given_file_text_lines = given_file_text_lines[1:]
    section = 'vertices'
    
    vertex_count = 0
    face_count = 0
    n_face_vertices = 0
    
    returning_indices_enabled = False
    
    vertex_count = result['n_vertices']
    if given_include_mesh:
        first_line_n_tokens = len(given_file_text_lines[0].split())
        result['vertices'] = xr.DataArray([tuple(map(float, x.split()[1:])) for x in given_file_text_lines[:vertex_count]], dims=('Vertex index', 'Dimension'), coords={'Vertex index': np.arange(vertex_count), 'Dimension': np.arange(first_line_n_tokens - 1)})
    given_file_text_lines = given_file_text_lines[vertex_count:]
    
    face_count = result['n_faces']
    if given_include_mesh:
        first_line_n_tokens = len(given_file_text_lines[0].split())
        result['faces'] = xr.DataArray([tuple(map(int, x.split()[3:])) for x in given_file_text_lines[:face_count]], dims=('Face index', 'Face vertex index'), coords={'Face index': np.arange(face_count), 'Face vertex index': np.arange(first_line_n_tokens - 3)})
        int_data = np.array([tuple(map(int, x.split()[:2])) for x in given_file_text_lines[:face_count]])
        if returning_indices_enabled:
            result['face_index_to_object_index'] = xr.DataArray(int_data[:, 1], dims=('Face index'), coords={'Face index': np.arange(face_count)})
            object_count = len(result['face_index_to_object_index'].unique())
            result['n_objects'] = object_count
            result['object_index_to_face_indices'] = {x: (result['face_index_to_object_index'] == x).find() for x in range(object_count)}
        result['faces_types'] = xr.DataArray([x.split()[2] for x in given_file_text_lines[:face_count]], dims=('Face index',), coords={'Face index': np.arange(face_count)})
    given_file_text_lines = given_file_text_lines[face_count:]
    
    result['face_data_first_line'] = list(map(int, given_file_text_lines[0].split()))
    result['n_face_data_variables'] = len(result['face_data_first_line']) - 1
    given_file_text_lines = given_file_text_lines[1:]
    
    result['variable_names'] = []
    result['face_data_variables_units'] = {}
    for i in range(result['n_face_data_variables']):
        given_file_text_line = given_file_text_lines[i]
        given_file_text_line_tokens = given_file_text_line.split(',')
        result['variable_names'].append(given_file_text_line_tokens[0])
        result['face_data_variables_units'][given_file_text_line_tokens[0]] = given_file_text_line_tokens[1]
    given_file_text_lines = given_file_text_lines[result['n_face_data_variables']:]
    
    first_line_n_tokens = len(given_file_text_lines[0].split())
    result['face_data'] = xr.DataArray([tuple(map(float, x.split()[1:])) for x in given_file_text_lines[:face_count]], dims=('Face index', 'Variable'), coords={'Face index': np.arange(face_count), 'Variable': result['variable_names']})
    
    return result


def parse_WSV_file(given_file_path):
    '''
    References:
    
    * 
    '''
    
    with open(given_file_path, 'r') as given_file:
        given_file_text_lines = given_file.readlines()
    
    column_names = None
    times = []
    values = []
    
    for given_file_text_line in given_file_text_lines:
        given_file_text_line = given_file_text_line.strip()
        if len(given_file_text_line) == 0:
            continue
        given_file_text_line_tokens = given_file_text_line.split()
        if given_file_text_line.startswith('#'):
            if column_names is None and len(times) == 0:
                column_names = given_file_text_line_tokens[1:]
            continue
        times.append(float(given_file_text_line_tokens[0]))
        values.append([float(x) for x in given_file_text_line_tokens[1:]])
    
    if column_names is None:
        result = pd.DataFrame(values, index=times, columns=column_names)
    else:
        result = pd.DataFrame(values, index=pd.Index(times, name=column_names[0]), columns=column_names[1:])
    
    return result


def parse_CSV_file(given_file_path):
    '''
    References:
    
    * 
    '''
    
    with open(given_file_path, 'r') as given_file:
        given_file_text_lines = given_file.readlines()
    
    column_names = None
    times = []
    values = []
    
    for given_file_text_line in given_file_text_lines:
        given_file_text_line = given_file_text_line.strip()
        if len(given_file_text_line) == 0:
            continue
        given_file_text_line_tokens = given_file_text_line.split(',')
        if column_names is None and len(times) == 0:
            if given_file_text_line.startswith('#'):
                column_names = given_file_text_line_tokens[1:]
            else:
                column_names = given_file_text_line_tokens
            continue
        times.append(float(given_file_text_line_tokens[0]))
        values.append([float(x) for x in given_file_text_line_tokens[1:]])
    
    if column_names is None:
        result = pd.DataFrame(values, index=times, columns=column_names)
    else:
        result = pd.DataFrame(values, index=pd.Index(times, name=column_names[0]), columns=column_names[1:])
    
    return result


vcell_column_name_pattern = r'\(Var=(.+)\)\s+\1'
vcell_column_name_pattern_program = re.compile(vcell_column_name_pattern)
vcell_column_name_pattern_program_matches = lambda x: vcell_column_name_pattern_program.fullmatch(x) is not None
vcell_column_name_simplify = lambda x: vcell_column_name_pattern_program.fullmatch(x).group(1)

def parse_SV_file(given_file_path):
    '''
    References:
    
    * 
    '''
    
    given_file_path_extension = os.path.splitext(given_file_path)[1][1:]
    if given_file_path_extension in ['wsv', 'txt', 'cdat', 'gdat']:
        result = parse_WSV_file(given_file_path)
    elif given_file_path_extension in ['csv']:
        result = parse_CSV_file(given_file_path)
        for i in range(len(result.columns)):
            result.columns.values[i] = vcell_column_name_simplify(result.columns.values[i])
    else:
        raise ValueError('Unknown file type with extension "{0}"'.format(given_file_path_extension))
    
    return result


def parse_DAT_file(given_file_path, given_include_mesh=False):
    '''
    References:
    
    * 
    '''
    
    geometry_name = given_file_path.split(os.path.sep)[-6]
    variable_name = os.path.splitext(os.path.splitext(given_file_path.split(os.path.sep)[-1])[0])[0]
    
    with open(given_file_path, 'r') as given_file:
        given_file_text_lines = given_file.readlines()
    
    times = []
    values = []
    
    for given_file_text_line in given_file_text_lines:
        given_file_text_line = given_file_text_line.strip()
        if len(given_file_text_line) == 0:
            continue
        if given_file_text_line.startswith('#'):
            continue
        given_file_text_line_tokens = given_file_text_line.split()
        times.append(float(given_file_text_line_tokens[0]))
        values.append(float(given_file_text_line_tokens[1]))
    
    result = pd.Series(values, index=pd.Index(times, name='Time'), name=(geometry_name, variable_name))
    
    return result


def parse_MCell_DAT_file(given_file_path, given_return_only_counts=False):
    '''
    References:
    
    * CellBlender Format: https://github.com/mcellteam/mcell/blob/5c36ac180b1f481b4ded30f4d2579e7fed8281e3/src/viz_output.c#L1317
    
    * https://github.com/mcellteam/mcell/blob/5c36ac180b1f481b4ded30f4d2579e7fed8281e3/src/viz_output.c#L2470
    
    * https://github.com/mcellteam/mcell/blob/5c36ac180b1f481b4ded30f4d2579e7fed8281e3/src/mcell_viz.c#L58
    * https://github.com/mcellteam/cellblender/blob/b2cd2b5fc277f6a26d0a53070784b61cce233fe8/engine_runner_combos/libMCell.cpp#L195
    * https://github.com/mcellteam/cellblender/blob/b2cd2b5fc277f6a26d0a53070784b61cce233fe8/engine_runner_combos/pure_python_sim.py
    '''
    
    debug_parse_MCell_DAT_file = False
    # debug_parse_MCell_DAT_file = True
    
    filename = os.path.split(given_file_path)[1]
    
    with open(given_file_path, 'rb') as given_file:
        mc_file_content = given_file.read()
    
    mc_big_endian = mc_file_content[0] == 0
    mc_byte_order = 'big' if mc_big_endian else 'little'
    mc_float_format = '>f' if mc_big_endian else '<f'
    def read_endian_bytes(given_offset, given_count=1):
        mc_result = mc_file_content[given_offset:given_offset + given_count]
        return mc_result
    
    def mc_to_int(given_bytes):
        return int.from_bytes(given_bytes, mc_byte_order)
    
    mc_file_content_length = len(mc_file_content)
    mc_file_content_offset = 0
    def mc_consume_bytes(given_count=1):
        nonlocal mc_file_content_offset
        mc_result = mc_file_content[mc_file_content_offset:mc_file_content_offset + given_count]
        mc_file_content_offset += given_count
        if debug_parse_MCell_DAT_file:
            print('mc_consume_bytes({}) returning {}'.format(given_count, mc_result))
        return mc_result
    
    def mc_ignore_bytes(given_count):
        nonlocal mc_file_content_offset
        mc_file_content_offset += given_count
        if debug_parse_MCell_DAT_file:
            print('mc_ignore_bytes({})'.format(given_count))
    
    def mc_consume_int(given_count=4):
        mc_result = int.from_bytes(mc_consume_bytes(given_count), mc_byte_order)
        if debug_parse_MCell_DAT_file:
            print('mc_consume_int({}) returning {}'.format(given_count, mc_result))
        return mc_result
    
    def mc_consume_float():
        mc_result = struct.unpack(mc_float_format, mc_consume_bytes(4))[0]
        if debug_parse_MCell_DAT_file:
            print('mc_consume_float() returning {}'.format(mc_result))
        return mc_result
    
    def mc_consume_floats(given_count):
        mc_result = np.array([mc_consume_float() for x in range(given_count)])
        if debug_parse_MCell_DAT_file:
            print('mc_consume_floats({}) returning {}'.format(given_count, mc_result))
        return mc_result
    
    def mc_consume_str():
        # return mc_consume_bytes(mc_consume_int(1).decode('ascii'))
        temp1 = mc_consume_int(1)
        temp2 = mc_consume_bytes(temp1)
        # No terminating null
        # mc_ignore_bytes(1)
        mc_result = temp2.decode('ascii')
        if debug_parse_MCell_DAT_file:
            print('mc_consume_str() with length {}, bytes {} returning {}'.format(temp1, temp2, mc_result))
        return mc_result
    
    mc_version = mc_consume_int()
    
    if mc_version != 1:
        raise NotImplementedError('CellBlender viz_output format %d not supported', mc_version)
    
    blocks_names = []
    blocks_molecule_types = []
    blocks_n_molecules = []
    if not given_return_only_counts:
        blocks_molecules_x = []
        blocks_molecules_y = []
        blocks_molecules_z = []
        blocks_molecules_orientation_x = []
        blocks_molecules_orientation_y = []
        blocks_molecules_orientation_z = []
    
    while mc_file_content_offset < mc_file_content_length:
        block_name = mc_consume_str()
        block_molecule_type = mc_consume_int(1)
        block_n_molecules = mc_consume_int() // 3
        
        blocks_names.append(block_name)
        blocks_molecule_types.append(block_molecule_type)
        blocks_n_molecules.append(block_n_molecules)
        
        if given_return_only_counts:
            mc_ignore_bytes(block_n_molecules * 3 * 4)
            if block_molecule_type == 1:
                mc_ignore_bytes(block_n_molecules * 3 * 4)
        else:
            block_molecules_xyz = mc_consume_floats(block_n_molecules * 3)
            blocks_molecules_x.append(block_molecules_xyz[0::3])
            blocks_molecules_y.append(block_molecules_xyz[1::3])
            blocks_molecules_z.append(block_molecules_xyz[2::3])
            if block_molecule_type == 1:
                # Surface molecules have orientation vectors
                block_molecules_orientation_xyz = mc_consume_floats(block_n_molecules * 3)
                blocks_molecules_orientation_z.append(block_molecules_orientation_xyz[2::3])
                blocks_molecules_orientation_x.append(block_molecules_orientation_xyz[0::3])
                blocks_molecules_orientation_y.append(block_molecules_orientation_xyz[1::3])
            else:
                blocks_molecules_orientation_x.append(None)
                blocks_molecules_orientation_y.append(None)
                blocks_molecules_orientation_z.append(None)
    
    result = {
        'filename': filename,
        'iteration': int(re.search('\.cellbin\.([0-9]+).dat', filename).group(1)),
        'variable_names': blocks_names,
        'n_molecules': pd.Series(blocks_n_molecules, index=blocks_names, name=filename),
        'molecule_type': pd.Series(blocks_molecule_types, index=blocks_names, name=filename),
        }
    if not given_return_only_counts:
        result['x'] = [pd.Series(v, name=u) for u, v in zip(blocks_names, blocks_molecules_x)]
        result['y'] = [pd.Series(v, name=u) for u, v in zip(blocks_names, blocks_molecules_y)]
        result['z'] = [pd.Series(v, name=u) for u, v in zip(blocks_names, blocks_molecules_z)]
        result['ox'] = {u: pd.Series(v, name=u) for u, v in zip(blocks_names, blocks_molecules_orientation_x) if v is not None}
        result['oy'] = {u: pd.Series(v, name=u) for u, v in zip(blocks_names, blocks_molecules_orientation_y) if v is not None}
        result['oz'] = {u: pd.Series(v, name=u) for u, v in zip(blocks_names, blocks_molecules_orientation_z) if v is not None}
        
        
        result['DataArray'] = {}
        for i, block_name in enumerate(blocks_names):
            """
            temp = {}
            temp['x'] = blocks_molecules_x[i]
            temp['y'] = blocks_molecules_y[i]
            temp['z'] = blocks_molecules_z[i]
            temp['ox'] = blocks_molecules_orientation_x[i]
            temp['oy'] = blocks_molecules_orientation_y[i]
            temp['oz'] = blocks_molecules_orientation_z[i]
            """
            block_n_molecules = blocks_n_molecules[i]
            temp = xr.DataArray(np.zeros((block_n_molecules, 6)) + np.nan, dims=('index', 'coord'), coords={'index': np.arange(block_n_molecules), 'coord': ['x', 'y', 'z', 'ox', 'oy', 'oz']})
            temp.loc[dict(coord='x')] = blocks_molecules_x[i]
            temp.loc[dict(coord='y')] = blocks_molecules_y[i]
            temp.loc[dict(coord='z')] = blocks_molecules_z[i]
            temp.loc[dict(coord='ox')] = blocks_molecules_orientation_x[i]
            temp.loc[dict(coord='oy')] = blocks_molecules_orientation_y[i]
            temp.loc[dict(coord='oz')] = blocks_molecules_orientation_z[i]
            result['DataArray'][block_name] = temp
    
    return result


def generate_and_simulate_generic_batch_id(reaction_network_name, downsampling, synthesis, simulation_end_time):
    '''
    Generate a data set ID used in `generate_and_simulate_generic_data_set_id` and `generate_and_simulate_generic_directory_name_pattern`.
    '''
    downsampling2 = float(downsampling)
    if synthesis == 'all':
        downsampling2 *= 1/4
    
    # return f"{reaction_network_name}_{downsampling2:.2e}_{simulation_end_time:.4e}sec"
    return f"{reaction_network_name}_{downsampling2:.2e}"


def generate_and_simulate_generic_data_set_id(reaction_network_name, downsampling, synthesis, simulation_end_time):
    '''
    Generate a data set ID used by `generate_and_simulate_analysis.py`.
    '''
    batch_id = generate_and_simulate_generic_batch_id(reaction_network_name, downsampling, synthesis, simulation_end_time)
    return f"{batch_id}_{simulation_end_time:.4e}sec"


def generate_and_simulate_generic_directory_name_pattern(reaction_network_name, downsampling, synthesis, simulation_end_time):
    '''
    Generate output directory name used by `generate_and_simulate_generic.m`.
    '''
    batch_id = generate_and_simulate_generic_batch_id(reaction_network_name, downsampling, synthesis, simulation_end_time)
    return f"img_{batch_id}_[0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9]_{simulation_end_time:.4e}sec"


def generate_and_simulate_generic_mcell_output_pattern(reaction_network_name, downsampling, synthesis, simulation_end_time):
    '''
    Generate output directory name used by `generate_and_simulate_generic.m`.
    '''
    return os.path.join(generate_and_simulate_generic_directory_name_pattern(reaction_network_name, downsampling, synthesis, simulation_end_time), 'cell[0-9]*', 'viz_data', 'seed_[0-9]*')


generate_and_simulate_generic_directory_name_geometry_seed_pattern = re.compile(r'img_(.*)_(' + float_pattern + r')_([0-9]+)_(' + float_pattern + r'sec)')
generate_and_simulate_generic_directory_name_program = re.compile(generate_and_simulate_generic_directory_name_geometry_seed_pattern)
def generate_and_simulate_generic_directory_name_geometry_seed(directory_name):
    '''
    Extract seed from data set ID used by `generate_and_simulate_analysis.py`.
    '''
    return generate_and_simulate_generic_directory_name_program.fullmatch(directory_name).group(3)


def parse_shape_space_coordinates_file(given_file_path, index_name='Shape space dimension', space_space_dimension_name_func=lambda x: 'Shape space dimension {0:d}'.format(x)):
    '''
    Parse output from CellOrganizer with `options.output.shape_space_coords = true`.
    '''
    
    filename = os.path.split(given_file_path)[1]
    geometry_name = given_file_path.split(os.path.sep)[-3]
    
    with open(given_file_path, 'r') as given_file:
        given_file_text = given_file.read()
    
    values = [float(x) for x in given_file_text.split()]
    
    space_space_dimension_names = [space_space_dimension_name_func(i) for i in range(len(values))]
    
    result = pd.Series(values, name=geometry_name, index=pd.Index(space_space_dimension_names, name=index_name))
    
    return result


def remove_vertices(given_vertices, given_faces, given_face_attributes, given_vertices_mask):
    given_faces_mask = given_vertices_mask[given_faces].all(axis=1)
    return remove_faces(given_vertices, given_faces, given_face_attributes, given_faces_mask)


def remove_faces(given_vertices, given_faces, given_face_attributes, given_faces_mask):
    # Remove faces
    given_faces = given_faces[given_faces_mask, :]
    for key, value in given_face_attributes.items():
        given_face_attributes[key] = given_face_attributes[key][given_faces_mask]
    # Remove unused vertices
    given_vertices_mask = np.zeros(given_vertices.shape[0], dtype='bool')
    given_vertices_mask[np.unique(given_faces)] = True
    new_vertex_indices = np.hstack(([0], np.cumsum(given_vertices_mask)[:-1]));
    given_vertices = given_vertices[given_vertices_mask, :]
    given_faces = new_vertex_indices[given_faces]
    return given_vertices, given_faces, given_face_attributes, given_vertices_mask


"""
def compute_UCD_statistics(given_all_UCD_info, given_save_mesh=False):
    '''
    '''
    
    given_UCD_info_keys = filename_sorted(given_all_UCD_info.keys())
    
    features = {}
    features['index'] = []
    pointwise_features = {}
    
    pointwise_feature_funcs = {}
    pointwise_feature_funcs['mean'] = np.mean
    pointwise_feature_funcs['max_minus_min'] = lambda x: np.max(x) - np.min(x)
    pointwise_feature_funcs['std'] = np.std
    pointwise_feature_funcs['std_over_mean'] = lambda x: np.std(x) - np.mean(x)
    pointwise_feature_funcs['median'] = np.median
    pointwise_feature_funcs['iqr'] = lambda x: np.quantile(x, 0.75) - np.quantile(x, 0.25)
    pointwise_feature_funcs_names = filename_sorted(pointwise_feature_funcs.keys())
    
    temporal_feature_funcs = {}
    temporal_feature_funcs['identity'] = lambda x: x
    temporal_feature_funcs['slope'] = lambda x: np.gradient(x, edge_order=2)
    def absolute_slope_falls_to_proportion(given_value, given_proportion):
        given_value_absolute_slope = np.absolute(np.gradient(given_value, edge_order=2))
        return given_value_absolute_slope <= given_value_absolute_slope.max() * given_proportion
    temporal_feature_funcs['abs_slope_falls_to_0.1_percent'] = lambda x: absolute_slope_falls_to_proportion(x, 1e-3)
    temporal_feature_funcs_names = filename_sorted(temporal_feature_funcs.keys())
    
    feature_combinations = list(product(pointwise_feature_funcs_names, temporal_feature_funcs_names))
    feature_combination_name_func = lambda x: ' '.join(x)
    # feature_combinations_map = {feature_combination_name_func(x): x for x in feature_combinations}
    
    for pointwise_feature_funcs_name in pointwise_feature_funcs_names:
        pointwise_features[pointwise_feature_funcs_name] = []
    
    for feature_combination in feature_combinations:
        feature_combination_name = feature_combination_name_func(feature_combination)
        features[feature_combination_name] = []
    
    
    # Process mesh
    
    given_UCD_info_key = given_UCD_info_keys[0]
    given_UCD_info = given_all_UCD_info[given_UCD_info_key]
    file_path = given_UCD_info['file_path']
    filename = given_UCD_info['filename']
    vertices = given_UCD_info['vertices']
    faces = given_UCD_info['faces']
    face_data = given_UCD_info['face_data']
    face_attributes = {'face_data': face_data}
    
    # Triangulate
    faces = np.vstack((faces[:, [0, 1, 2]], faces[:, [2, 3, 0]]))
    for key, value in face_attributes.items():
        face_attributes[key] = np.vstack((value,) * 2)
    
    # Remove faces where the concentration is zero (not participating in reaction or diffusion)
    faces_mask = (face_attributes['face_data'] != 0).ravel()
    vertices, faces, face_attributes, vertices_mask = remove_faces(vertices, faces, face_attributes, faces_mask)
    
    # Create Trimesh and save
    mesh = trimesh.Trimesh(vertices, faces, face_attributes=face_attributes, process=False)
    if given_save_mesh:
        # export_file_path = os.path.join(os.path.split(file_path)[0], filename + '.obj')
        export_file_path = os.path.join(base_filename, filename + '.obj')
        mesh.export(export_file_path)
    
    print('vertices.shape {}'.format(vertices.shape))
    print('faces.shape {}'.format(faces.shape))
    print('face_data.shape {}'.format(face_data.shape))
    
    
    # Compute pointwise features
    
    for i, given_UCD_info_key in enumerate(given_UCD_info_keys):
        given_UCD_info = given_all_UCD_info[given_UCD_info_key]
        
        face_data = given_UCD_info['face_data']
        
        # # First principal component of vertex coordinates (PC1) as substitute for principal axis of cell shape (would require watertight mesh)
        # vertices_pca = sklearn.decomposition.PCA(n_components=3)
        # vertices_pca.fit(vertices)
        
        features['index'].append(i)
        
        for pointwise_feature_funcs_name in pointwise_feature_funcs_names:
            pointwise_features[pointwise_feature_funcs_name].append(pointwise_feature_funcs[pointwise_feature_funcs_name](face_data))
    
    features['index'] = np.array(features['index'])
    
    # Compute temporal features
    
    for feature_combination in feature_combinations:
        feature_combination_name = feature_combination_name_func(feature_combination)
        pointwise_feature_funcs_name = feature_combination[0]
        temporal_feature_funcs_name = feature_combination[1]
        pointwise_feature = pointwise_features[pointwise_feature_funcs_name]
        features[feature_combination_name] = temporal_feature_funcs[temporal_feature_funcs_name](pointwise_feature)
    
    return features
"""




def load_shape_space_coordinates(coords_base_dir, coords_dir_pattern=None):
    coords_base_dir = os.path.expanduser(coords_base_dir)
    if coords_dir_pattern is not None:
        coords_base_dir = os.path.join(coords_base_dir, coords_dir_pattern)
    # print('coords_base_dir', coords_base_dir)
    coords_dir_list = glob.glob(os.path.join(coords_base_dir))
    coords_dir_list = list(coords_dir_list)
    coords_dir_list = map(lambda x: os.path.join(x, 'cell1') if not x.endswith('cell1') else x, coords_dir_list)
    coords_dir_list = filename_sorted(coords_dir_list)
    coords_dir_list = list(coords_dir_list)

    all_shape_space_coordinates = []

    for coords_dir in coords_dir_list:
        # Parse files
        
        coords_path = os.path.join(coords_dir, 'shape_space_coords.txt');
        # print('coords_path')
        # print(coords_path)
        shape_space_coordinates = parse_shape_space_coordinates_file(coords_path)
        all_shape_space_coordinates.append(shape_space_coordinates)
        
    all_shape_space_coordinates = pd.concat(all_shape_space_coordinates, axis=1)
    all_shape_space_coordinates.columns.names = ['Geometry']

    # print('all_shape_space_coordinates')
    # print(all_shape_space_coordinates)
    
    return all_shape_space_coordinates



def plot_and_save_mean_concentrations(given_mean_concentrations, given_prefix):
    given_geometry_names = given_mean_concentrations.coords['Geometry name'].values
    given_variable_names = given_mean_concentrations.coords['Variable'].values
    
    if should_plot_pair_mean_concentrations:
        for given_geometry_name, given_variable_name in product(given_geometry_names, given_variable_names):
            given_geometry_name_filename = os.path.join(base_filename, 'mean_concentrations {0} {1} {2}.png'.format(given_prefix, given_geometry_name, given_variable_name))
            if not os.path.exists(given_geometry_name_filename):
                plt.close('all')
                plt.figure(figsize=(figure_width, figure_height))
                plt.plot(given_mean_concentrations.loc[{'Geometry name': given_geometry_name, 'Variable': given_variable_name}])
                plt.savefig(given_geometry_name_filename)
    
    if should_plot_simulation_mean_concentrations:
        # Plot given_mean_concentrations for all variables in each simulation
        
        loaded_geometry_names = mean_concentrations.coords['Geometry name'].values
        
        for given_geometry_name in given_geometry_names:
            simulation_data = given_mean_concentrations.loc[{'Geometry name': given_geometry_name}].reset_coords(drop=True)
            simulation_data = simulation_data.rename('{0}'.format(given_geometry_name))
            plot_xarray(simulation_data, ['Variable'], 'concentrations in simulation', 'Concentration')
    
    if should_plot_variable_mean_concentrations:
        # Plot given_mean_concentrations for each variable in all simulations
        
        for variable_name in given_variable_names:
            variable_data = given_mean_concentrations.loc[{'Variable': variable_name}].reset_coords(drop=True)
            variable_data = variable_data.rename('{0}'.format(variable_name))
            plot_xarray(variable_data, ['Geometry name'], 'concentrations across simulations for species', 'Concentration', given_use_legend=False)


# def plot_and_save(value, variable_name_tuples, prefix):
def plot_and_save_xarray_combinations(data, x_variable_name, distinct_plot_variable_names, base_filename, figure_size=(4, 3), should_overwrite=False):
    '''
    Plot values vs x_variable_name for combinations of categorical variables.
    
    distinct_plot_variable_names: Create a separate plot for combinations of these variables
    '''
    x_variable_coords = data.coords[x_variable_name].values
    
    # data_other_variables_names = list(filter(lambda x: x not in (x_variable_name, y_variable_name), data.dims))
    # data_other_variables_coords = [data.coords[x].values for x in data_other_variables_names]
    distinct_variables_coords = [data.coords[x].values for x in distinct_plot_variable_names]
    for distinct_variables_coord in product(*distinct_variables_coords):
        distinct_variables_coords_str = ','.join([str(x) for x in distinct_variables_coord])
        geometry_name_filename = '{} {} vs {} for {}.png'.format(base_filename, data.name, x_variable_name, distinct_variables_coords_str)
        if not os.path.exists(geometry_name_filename) or should_overwrite:
            plt.close('all')
            plt.figure(figsize=figure_size)
            distinct_variables_coord_dict = {x: y for x, y in zip(distinct_plot_variable_names, distinct_variables_coord)}
            temp = data.loc[distinct_variables_coord_dict].reset_coords(drop=True)
            plt.plot(temp)
            plt.title(distinct_variables_coords_str)
            plt.xlabel(x_variable_name)
            plt.ylabel(data.name)
            # plot_xarray(temp, ['Geometry name'], '?', 'Concentration', use_legend=False)
            plt.savefig(geometry_name_filename)
    
    """
    for given_variable_name_tuple in given_variable_name_tuples:
        for given_xarray_tuple in product(*given_variable_name_tuple):
            given_geometry_name_filename = os.path.join(base_filename, 'mean_concentrations {0} {1}.png'.format(given_prefix, ' '.join([str(x) for x in given_xarray_tuple])))
            if not os.path.exists(given_geometry_name_filename):
                # plt.close('all')
                # plt.figure(figsize=(figure_width, figure_height))
                # plt.plot(given_xarray.loc[{x: y for x, y in zip(given_variable_name_tuple, given_xarray_tuple)}])
                # plot_xarray(given_xarray.loc[{x: y for x, y in zip(given_variable_name_tuple, given_value_tuple)}, ['Geometry name'], '?', 'Concentration', given_use_legend=False)
                plt.savefig(given_geometry_name_filename)
    
    if should_plot_simulation_mean_concentrations:
        # Plot given_mean_concentrations for all variables in each simulation
        
        loaded_geometry_names = mean_concentrations.coords['Geometry name'].values
        
        for given_geometry_name in given_geometry_names:
            simulation_data = given_mean_concentrations.loc[{'Geometry name': given_geometry_name}].reset_coords(drop=True)
            simulation_data = simulation_data.rename('{0}'.format(given_geometry_name))
            plot_xarray(simulation_data, ['Variable'], 'concentrations in simulation', 'Concentration')
    
    if should_plot_variable_mean_concentrations:
        # Plot given_mean_concentrations for each variable in all simulations
        
        for variable_name in given_variable_names:
            variable_data = given_mean_concentrations.loc[{'Variable': variable_name}].reset_coords(drop=True)
            variable_data = variable_data.rename('{0}'.format(variable_name))
            plot_xarray(variable_data, ['Geometry name'], 'concentrations across simulations for species', 'Concentration', given_use_legend=False)
    """



def process_cellorganizer_mcell_output(data_base_dir):
    data_base_dir = os.path.expanduser(data_base_dir)
    print('data_base_dir', data_base_dir)
    data_dir_list = glob.glob(os.path.join(data_base_dir, 'img_LotkaVolterra_[0-9]*[0-9]'))
    coords_dir_list = list(data_dir_list)
    print('len(data_dir_list)', len(data_dir_list))
    # print('\nDEBUG\n'); data_dir_list = data_dir_list[:5]
    data_dir_list = map(lambda x: os.path.join(x, 'cell1', 'react_data', 'cell'), data_dir_list)
    data_dir_list = map(lambda x: glob.glob(os.path.join(x, 'seed_[0-9]*[0-9]')), data_dir_list)
    data_dir_list = map(lambda x: x[0] if len(x) > 0 else None, data_dir_list)
    data_dir_list = filter(lambda x: x is not None, data_dir_list)
    data_dir_list = filename_sorted(data_dir_list)
    data_dir_list = list(data_dir_list)
    # print('data_dir_list', data_dir_list)

    all_reaction_data = []
    all_shape_space_coordinates = []
    all_data_dir_names = []

    for data_dir in data_dir_list:
        # Parse files
        
        file_path_list = glob.glob(os.path.join(data_dir, '*.dat'))
        if len(file_path_list) == 0:
            continue
        data_dir_name = data_dir.split(os.path.sep)[-5]
        file_path_list = filename_sorted(file_path_list)
        dir_DAT_info = []
        for i, file_path in enumerate(file_path_list):
            DAT_info = parse_DAT_file(file_path, i == 0)
            dir_DAT_info.append(DAT_info)
        
        all_data_dir_names.append(data_dir_name)
        
        # data_dir_data = pd.merge(dir_DAT_info)
        data_dir_data = pd.concat(dir_DAT_info, axis=1)
        data_dir_data.columns.names = ['Geometry', 'Reaction data']
        all_reaction_data.append(data_dir_data)
        
        # coords_dir = os.path.join(data_dir.split(os.path.sep)[:-3])
        coords_dir = os.path.join(data_dir, '..', '..', '..')
        coords_path = os.path.join(coords_dir, 'shape_space_coords.txt');
        # print('coords_path')
        # print(coords_path)
        shape_space_coordinates = parse_shape_space_coordinates_file(coords_path)
        # print('shape_space_coordinates')
        # print(shape_space_coordinates)
        all_shape_space_coordinates.append(shape_space_coordinates)
        # raise
        

    def get_level_values_list(given_multiindex, given_level):
        given_multiindex_keys = [x.columns.get_level_values(given_level) for x in given_multiindex]
        given_multiindex_keys2 = pd.Index([])
        for x in given_multiindex_keys:
            given_multiindex_keys2 |= x
        given_multiindex_keys = given_multiindex_keys2
        given_multiindex_keys = given_multiindex_keys.unique()
        return given_multiindex_keys

    simulation_keys = get_level_values_list(all_reaction_data, 0)
    reaction_data_keys = get_level_values_list(all_reaction_data, 1)
    reaction_data_species_keys = list(filter(lambda x: re.fullmatch('s[0-9]+', x), reaction_data_keys))
    reaction_data_reaction_keys = list(filter(lambda x: re.fullmatch('r[0-9]+', x), reaction_data_keys))
    # all_reaction_data = pd.merge(all_reaction_data)
    # raise NotImplementedError
    all_reaction_data = pd.concat(all_reaction_data, axis=1)
    # all_reaction_data = pd.concat(all_reaction_data, axis=1, names=['Reaction output', 'Geometry'], keys=[reaction_data_keys, simulation_keys])
    print('all_reaction_data.columns.names'); print(all_reaction_data.columns.names)
    # Make reaction data keys level 0
    all_reaction_data.columns = all_reaction_data.columns.swaplevel(0, 1)
    # reaction_data_keys = all_reaction_data.columns.get_level_values(0).unique()
    # simulation_keys = all_reaction_data.columns.get_level_values(1).unique()


    # print('all_reaction_data')
    # print(all_reaction_data)

    # all_shape_space_coordinates = pd.merge(all_shape_space_coordinates)
    # raise(NotImplementedError)  
    # all_shape_space_coordinates = pd.concat(all_shape_space_coordinates, axis=1, names=['Geometry'])
    all_shape_space_coordinates = pd.concat(all_shape_space_coordinates, axis=1)
    all_shape_space_coordinates.columns.names = ['Geometry']

    # print('all_shape_space_coordinates')
    # print(all_shape_space_coordinates)

    all_reaction_data_values = all_reaction_data.values
    n_times = all_reaction_data_values.shape[0]
    window = sp.signal.blackman(n_times)
    all_reaction_data_fourier = sp.fftpack.fft(all_reaction_data.values * window[:, None])
    all_reaction_data_fourier = np.abs(all_reaction_data_fourier)
    # for x in all_reaction_data:
        # print('x')
        # print(x)
        # raise NotImplementedError



def analyze_cellorganizer_mcell_output():
    # Plot for each of the reaction data

    regression_random_state = 106435
    model_steps = []
    model_steps.append(('standardization', StandardScaler(with_mean=True, with_std=True)))
    # model_steps.append(('regression', LinearRegression(fit_intercept=True, normalize=False)))
    model_steps.append(('regression', TheilSenRegressor(fit_intercept=True, max_subpopulation=100, random_state=regression_random_state)))
    model = Pipeline(model_steps)
    model_type = ','.join([x[0] for x in model.steps])
    n_splits = 5
    n_repeats = 10
    model_cv_random_state = 349457
    model_cv = RepeatedKFold(n_splits=n_splits, n_repeats=n_repeats, random_state=model_cv_random_state)

    # model_scoring = 'neg_mean_squared_error'
    model_scoring = 'r2'


    print('model'); print(model)
    print('model_cv'); print(model_cv)
    print('model_scoring {0}'.format(model_scoring))


    # default_font_size = 10
    # default_font_size = 16
    default_font_size = 36
    # default_marker_size = 6
    default_marker_size = 12
    matplotlib.rc('font', size=default_font_size)
    matplotlib.rc('lines', markersize=default_marker_size)
    matplotlib.rcParamsDefault['lines.markersize']

    # Y_trim = 1e-2
    Y_trim = 2.5e-2
    # Y_trim = 5e-2
    # Y_trim = 1e-1

    print('Y_trim {0}'.format(Y_trim))


    shape_space_dimension_prefix = 'Shape space dimension'
    shape_space_dimension_prefix2 = 'Dim.'
    rms_normalized_derivative_suffix = 'RMS derivative over initial count'
    rms_normalized_derivative_suffix2 = 'norm. RMS deriv.'
    def name_display_transform_func(x):
        if x.startswith(shape_space_dimension_prefix):
            # x = os.linesep.join(['Dim', x.split(shape_space_dimension_prefix)[1]])
            # x = '{0} {1}'.format(shape_space_dimension_prefix2, x.split(shape_space_dimension_prefix)[1].trim())
            x = x.replace(shape_space_dimension_prefix, shape_space_dimension_prefix2)
        elif x.endswith(rms_normalized_derivative_suffix):
            x = x.replace(rms_normalized_derivative_suffix, rms_normalized_derivative_suffix2)
        return x


    # raise NotImplementedError('Compute regressions, correct with sm.stats.multipletests, plot grid directly through pyplot')

    def annotation_func(x, y, **kwargs):
        if model_type == 'OLS':
            model_results = sm.OLS(y.values, sm.add_constant(x.values)).fit()
            plt.plot(x, y)
            axes = plt.gca()
            print('model_results.params'); print(model_results.params)
            sm.graphics.abline_plot(model_results=model_results)
            f_test_results = model_results.f_test(np.identity(2))
            print('f_test_results'); print(f_test_results)
            print('dir(f_test_results)'); print(dir(f_test_results))
            axes.annotate('F-test p = {0:.3e}'.format(model_results.f_test(np.identity(2)).pvalue))
        elif model_type == 'OLS MV':
            # pair_pipeline = sklearn.pipeline.make_pipeline(*pipeline_steps)
            # pair_pipeline.fit(x, y)
            raise NotImplementedError

    # color_arg_min = (1, 1, 1)
    # color_arg_max = (0, 0, 0)
    def scatter_with_continuous_hue_func(x, y, **kwargs):
        # https://stackoverflow.com/questions/44641669/scatterplot-with-point-colors-representing-a-continuous-variable-in-seaborn-face
        color_arg = kwargs['color']
        kwargs.pop('color')
        # print('color_arg {0}'.format(color_arg))
        global color_arg_min
        global color_arg_max
        color_arg_min = np.min([color_arg_min, color_arg], axis=0)
        color_arg_max = np.max([color_arg_max, color_arg], axis=0)
        plt.scatter(x, y, c=color_arg, **kwargs)

    # palette = sns.cubehelix_palette(n_colors=1000, start=1.25, rot=0.4, reverse=False)
    # palette = sns.cubehelix_palette(n_colors=1000, start=1.25, rot=0.75, dark=0.5, light=0.5, reverse=False)
    # palette = sns.cubehelix_palette(n_colors=1000, start=1.25, rot=270, reverse=False)
    palette = sns.mpl_palette('viridis', n_colors=1000)
    cmap = matplotlib.colors.ListedColormap(palette)

    all_reaction_data_times = all_reaction_data.index.to_series()
    all_reaction_data_time_step = ((all_reaction_data_times.shift(-1) - all_reaction_data_times.shift(1)) / 2).median()
    # for reaction_data_key in reaction_data_keys:
    for species_key in reaction_data_species_keys:
        species_data = all_reaction_data[species_key]
        species_data_central_difference_dataframe = (species_data.shift(-1) - species_data.shift(1)) / 2
        species_data_derivative_dataframe = species_data_central_difference_dataframe / all_reaction_data_time_step
        species_data_initial_count = species_data.iloc[0]

        # RMS reaction rate (central difference divided by time step) divided by initial molecule count and simulation time
        species_data_derivative_normalized_rms = ((species_data_derivative_dataframe / species_data_initial_count)**2).mean(axis=0).pow(1/2)
        species_data_derivative_normalized_rms_name = '{0} RMS derivative over initial count'.format(species_key)
        species_data_derivative_normalized_rms.name = species_data_derivative_normalized_rms_name
        species_data_derivative_normalized_rms_dataframe = pd.DataFrame(species_data_derivative_normalized_rms)

        # plot_data = pd.merge(all_shape_space_coordinates, species_data_derivative_normalized_rms_dataframe.T, on='Geometry')
        plot_data = pd.concat([all_shape_space_coordinates.T, species_data_derivative_normalized_rms_dataframe], axis=1)
        plot_data.dropna(inplace=True)
        
        # print('plot_data.shape {0}'.format(plot_data.shape)); raise # Debug
        
        
        # Regression of normalized species derivative onto shape space coordinates
        
        X_columns = list(filter(lambda x: x != species_data_derivative_normalized_rms_name, plot_data.columns))
        Y_column = species_data_derivative_normalized_rms_name
        X = plot_data[X_columns]
        Y = plot_data[Y_column]
        model_cv_scores = cross_val_score(model, X, Y, cv=model_cv, scoring=model_scoring)
        scores_mean = model_cv_scores.mean()
        scores_std = model_cv_scores.std()
        print('\'{0}\' vs. \'{1}\''.format(species_data_derivative_normalized_rms_name, 'all_shape_space_coordinates.T'))
        # print('    scores_mean {1:>12.4e}, scores_std {2:>12.4e}'.format(species_key, scores_mean, scores_std))
        # print('    scores_mean {1:>12.4f}, scores_std {2:>12.4e}'.format(species_key, scores_mean, scores_std))
        # print('    scores_mean {1:>12.4f}, scores_std {2:>12.4f}'.format(species_key, scores_mean, scores_std))
        print('    scores_mean {1:>7.4f}, scores_std {2:>7.4f}'.format(species_key, scores_mean, scores_std))
        
        Y_limits = (Y.quantile(Y_trim), Y.quantile(1 - Y_trim))
        
        palette2 = palette
        # palette2 = sns.mpl_palette('viridis', n_colors=len(Y))
        cmap2 = matplotlib.colors.ListedColormap(palette2)
        
        
        """
        # Scatter and KDE plots
        
        g = sns.PairGrid(X, height=3, aspect=1.6)
        g.map_upper(plt.scatter)
        g.map_lower(sns.kdeplot)
        g.map_diag(sns.kdeplot)
        for extension in ['png', 'svg']:
            figure_path = os.path.join(base_filename, '{0}_{1}_pairgrid.{2}'.format(model_type, species_key, extension))
            g.savefig(figure_path)
        del g.fig
        del g
        plt.close('all')
        """
        
        
        
        # Scatter plots with coloring
        
        n_X_columns = len(X.columns)
        height = 3
        aspect = 1.6
        fig = plt.figure(figsize=(height * aspect * (n_X_columns - 1), height * (n_X_columns - 1)))
        axis_grid = {x: {x: None for x in X.columns} for x in X.columns}
        for i1, c1 in enumerate(X.columns):
            for i2, c2 in enumerate(X.columns):
                if i2 >= i1:
                    continue
                fig_subplot = plt.subplot(n_X_columns - 1, n_X_columns - 1, (i1 - 1) * (n_X_columns - 1) + i2 + 1)
                axis_grid[c2][c1] = fig_subplot.scatter(X[c1], X[c2], c=Y, vmin=Y_limits[0], vmax=Y_limits[1], cmap=cmap2)
                fig_subplot.set_frame_on(False)
                fig_subplot.set_xticks([])
                fig_subplot.set_yticks([])
                if i1 == n_X_columns - 1:
                    plt.xlabel(name_display_transform_func(c2))
                if i2 == 0:
                    plt.ylabel(name_display_transform_func(c1))
        
        # https://stackoverflow.com/questions/44641669/scatterplot-with-point-colors-representing-a-continuous-variable-in-seaborn-face
        fig.subplots_adjust(right=0.92)
        # colorbar_axes = fig.add_axes([0.94, 0.25, 0.02, 0.6])
        colorbar_axes = fig.add_axes([0.94, 0.2, 0.02, 0.6])
        # colorbar_axes.set_ylabel(name_display_transform_func(Y.name))
        # colorbar_axes.set_title(name_display_transform_func(Y.name))
        # fig.suptitle(name_display_transform_func(Y.name))
        fig.suptitle(Y.name)
        # colorbar_points = plt.scatter([], [], c=[], vmin=Y.min(), vmax=Y.max(), cmap=cmap2)
        colorbar_points = plt.scatter([], [], c=[], vmin=Y_limits[0], vmax=Y_limits[1], cmap=cmap2)
        fig.colorbar(colorbar_points, cax=colorbar_axes)
        for extension in ['png', 'svg']:
            figure_path = os.path.join(base_filename, '{0}_{1}_pairgrid_normalized_derivative.{2}'.format(model_type, species_key, extension))
            plt.savefig(figure_path)
        del fig, fig_subplot, axis_grid, colorbar_axes, colorbar_points
        plt.close('all')
        
        
        
        """
        # Regression with individual coordinates
        
        # g = sns.PairGrid(plot_data, x_vars=all_shape_space_coordinates.index, y_vars=[species_data_derivative_normalized_rms.name], height=3, aspect=1.6)
        g = sns.PairGrid(plot_data, x_vars=X.columns, y_vars=Y.name, height=3, aspect=1.6)
        g.map(sns.regplot)
        trim = 1e-2
        g.set(ylim=species_data_derivative_normalized_rms.quantile([trim, 1 - trim]))
        for extension in ['png', 'svg']:
            figure_path = os.path.join(base_filename, '{0}_{1}_individual_regression.{2}'.format(model_type, species_key, extension))
            g.savefig(figure_path)
        del g.fig
        del g
        plt.close('all')
        """





