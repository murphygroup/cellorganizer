'''
Utilities for xarray.

Taraz Buck
2021-10-24
'''


from enable_interactive_debugger import sys


# import os
# import glob
# import fnmatch
# from math import *

# import numpy as np
import xarray as xr



xarray_forbidden_characters = []
xarray_forbidden_characters.append(' ')
xarray_forbidden_characters_encode_func = lambda x, y: x.replace(y, '[ord:{0:d}]'.format(ord(y)))
xarray_forbidden_characters_decode_func = lambda x, y: x.replace('[ord:{0:d}]'.format(ord(y)), y)

def reformat_xarray_for_saving(given_value):
    return_value = given_value.copy()
    for coord_name, coord in return_value.coords.items():
        coord_name_reformatted = str(coord_name)
        for xarray_forbidden_character in xarray_forbidden_characters:
            coord_name_reformatted = xarray_forbidden_characters_encode_func(coord_name_reformatted, xarray_forbidden_character)
        return_value = return_value.rename({coord_name: coord_name_reformatted})
    return return_value

def reformat_xarray_from_saved(given_value):
    return_value = given_value.copy()
    for coord_name_reformatted, coord in return_value.coords.items():
        coord_name = str(coord_name_reformatted)
        for xarray_forbidden_character in xarray_forbidden_characters:
            coord_name = xarray_forbidden_characters_decode_func(coord_name, xarray_forbidden_character)
        return_value = return_value.rename({coord_name_reformatted: coord_name})
    return return_value

def get_xarray_filename(given_base_filename, given_key):
    return '{0} {1}.nc'.format(given_base_filename, given_key)

def save_xarray(given_base_filename, **kwargs):
    for given_key, given_value in kwargs.items():
        given_value = reformat_xarray_for_saving(given_value)
        given_value.to_netcdf(get_xarray_filename(given_base_filename, given_key))

def load_xarray(given_base_filename, given_keys):
    return_value = {}
    for given_key in given_keys:
        return_value[given_key] = xr.open_dataarray(get_xarray_filename(given_base_filename, given_key))
        return_value[given_key] = reformat_xarray_from_saved(return_value[given_key])
    return return_value


def plot_xarray(given_xarray, given_lines_dimensions_names, given_filename_prefix, given_title=None, given_y_axis_label=None, given_use_legend=True, given_figure_size=(4, 3)):
    import matplotlib
    matplotlib.use('agg')
    import matplotlib.pyplot as plt
    import seaborn as sns
    
    # print('given_xarray'); print(given_xarray)
    # print('given_lines_dimensions_names'); print(given_lines_dimensions_names)
    # print('given_filename_prefix'); print(given_filename_prefix)
    # print('given_title'); print(given_title)
    # print('given_y_axis_label'); print(given_y_axis_label)
    # print('given_use_legend'); print(given_use_legend)
    
    given_lines_coords = [given_xarray.coords[x].values for x in given_lines_dimensions_names]
    given_lines_coords = product(*given_lines_coords)
    given_lines_coords = list(given_lines_coords)
    given_lines_coords_dicts = [dict(zip(given_lines_dimensions_names, x)) for x in given_lines_coords]
    given_lines_names = [','.join(str(x)) for x in given_lines_coords]
    
    # palette = sns.mpl_palette('viridis', n_colors=n_variables)
    # palette = distinctipy.get_colors(n_variables, exclude_colors=[(1, 1, 1)], pastel_factor=0.2, n_attempts=10000, colorblind_type='Deuteranomaly')
    # palette = color_tol.qualitative(n_variables)

    # cmap = matplotlib.colors.ListedColormap(palette)
    
    fig = plt.figure(figsize=given_figure_size)
    
    for given_line_name, given_line_coords_dict in zip(given_lines_names, given_lines_coords_dicts):
        given_line_data = given_xarray.loc[given_line_coords_dict].reset_coords(drop=True)
        
        # print(given_line_data); raise NotImplementedError
        
        # plt.plot(given_line_data, '-', color=palette[i])
        plt.plot(given_line_data)
        # raise NotImplementedError
    
    plt.xlabel(given_line_data.coords[given_line_data.dims[0]].name)
    if given_y_axis_label is None:
        given_y_axis_label = given_filename_prefix
    plt.ylabel(given_y_axis_label)
    if given_use_legend:
        # plt.legend(given_lines_names)
        # plt.legend(given_lines_names, 'right')
        # plt.legend(given_lines_names, 'right', ncol=3)
        plt.legend(given_lines_names, ncol=3)

    figure_filename = '{0} {1}'.format(given_filename_prefix, given_xarray.name)
    # figure_title = given_xarray.name
    # figure_title = figure_filename
    figure_title = given_title
    if given_title is None:
        figure_title = ''
    # fig.suptitle(figure_title)
    plt.title(figure_title)
    for extension in ['png']:
    # for extension in ['png', 'svg']:
        # figure_path = os.path.join(base_filename, '{0}.{1}'.format(figure_filename, extension))
        figure_path = os.path.join('{0}.{1}'.format(figure_filename, extension))
        plt.savefig(figure_path)
    del fig
    plt.close('all')


