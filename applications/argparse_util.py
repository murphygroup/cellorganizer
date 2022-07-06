'''
Utilities for argparse.

Taraz Buck
2021-10-12
'''


from enable_interactive_debugger import sys


import os
import glob
# import fnmatch
# from math import *
import argparse

# import numpy as np
# import xarray as xr



# import socket
# import tempfile
# from datetime import datetime
# import time

def remove_prepended_arguments(args):
    '''
    IPython sometimes returns the argv it was called with instead of one beginning with this script's name
    '''
    for i, arg in enumerate(args):
        if arg == __file__:
            return args[i:]
    return args

def nonnegative_int(x):
    x = int(x)
    if x < 0:
        raise argparse.ArgumentTypeError('Argument must be a nonnegative integer')
    return x

def positive_int(x):
    x = int(x)
    if x <= 0:
        raise argparse.ArgumentTypeError('Argument must be a positive integer')
    return x

def positive_float(x):
    x = float(x)
    if x <= 0:
        raise argparse.ArgumentTypeError('Argument must be a positive float')
    return x

def bool_strict(x):
    if x.lower() == 'true':
        return True
    elif x.lower() == 'false':
        return False
    else:
        raise argparse.ArgumentTypeError('Argument must be true or false')

# memory_spec_program = re.compile(r'[0-9]+[KMG]')

def existing_path(x):
    x = str(x)
    if not os.path.exists(x):
        print('type(x)', type(x))
        raise argparse.ArgumentTypeError(f'Argument must be an existing path (given {repr(x)})')
    return x

def existing_glob(x):
    x = str(x)
    if len(glob.glob(x)) == 0:
        print('type(x)', type(x))
        raise argparse.ArgumentTypeError(f'Argument must be an existing path or glob pattern matching at least one existing path (given {repr(x)})')
    return x

def list_with_type_func(x, type_func):
    x2 = []
    for y in x:
        x2.append(type_func(y))
    return x2

def list_with_type_func_func(type_func):
    return lambda x: list_with_type_func(x, type_func=type_func)

def single_or_list_with_type_func(x, type_func):
    try:
        return [type_func(x)]
    except argparse.ArgumentError:
        pass
    x2 = []
    for y in x:
        x2.append(type_func(y))
    return x2

def single_or_list_with_type_func_func(type_func):
    return lambda x: single_or_list_with_type_func(x, type_func=type_func)

# list_of_existing_paths = list_with_type_func_func(existing_path)
"""
def list_of_existing_paths(x):
    if isinstance(x, str):
        x = [x]
    x2 = []
    for y in x:
        x2.append(existing_path(y))
    return x2

def list_of_positive_integers(x):
    if isinstance(x, str):
        x = [x]
    x2 = []
    for y in x:
        x2.append(existing_path(y))
    return x2
"""
list_of_existing_paths = single_or_list_with_type_func_func(existing_path)
list_of_positive_integers = single_or_list_with_type_func_func(positive_int)
