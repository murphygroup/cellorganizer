'''
Utilities for parallel operations on a cluster.

Taraz Buck
2021-10-10
'''


from enable_interactive_debugger import sys


import os
import glob
# import fnmatch
from math import *
import types
import traceback
import inspect
import collections

import numpy as np
import xarray as xr

import xarray_util as xru



# import socket
import tempfile
from datetime import datetime
import time

def chunk_start(final_name, should_block=False, max_block_time=float('inf')):
    '''
    References:
    
    * http://www.dwheeler.com/secure-programs/Secure-Programs-HOWTO/avoid-race.html
    * https://docs.python.org/3/library/tempfile.html#tempfile.TemporaryFile
    * https://github.com/python/cpython/blob/975d10a4f8f5d99b01d02fc5f99305a86827f28e/Modules/_io/fileio.c#L305
    * https://pypi.org/project/atomicwrites/1.4.0/
    '''
    
    # final_root, final_extension = os.path.splitext(final_name)
    # temp_extension = '.tmp' if final_extension != '.tmp' else '.tmp.tmp'
    # if final_extension == '.tmp':
    # final_extension = fname_ext
    # temp_name = final_root + '.tmp'
    
    can_start = False
    final_exists = False
    temp_name = final_name + '.tmp'
    # temp_handle = None
    
    return_value_func = lambda: (can_start, final_exists, temp_name)
    if not os.path.exists(temp_name) and os.path.exists(final_name):
        final_exists = True
        return return_value_func()
    
    """
    # Create unique name using job ID and hostname
    hostname = socket.gethostname()
    job_id = os.environ.get('SLURM_JOB_ID', str(os.getpid()))
    unique_name = f'{final_name}_{hostname}_{job_id}.tmp'
    """
    
    """
    # Create unique file without race conditions
    final_name_head, final_name_tail = os.path.split(final_name)
    hard_link_count_func = lambda x: os.stat(x, follow_symlinks=False)
    unique_name = tempfile.mkstemp(suffix='.tmp', dir=final_name_head)
    hard_link_count_before = hard_link_count_func(unique_name)
    """
    
    # Attempt to create the temp file. If successful, return with can_start = True. Otherwise, wait or return with can_start = False.
    block_start_time = datetime.now()
    still_blocking = True
    while still_blocking:
        # Python open(file, 'x') calls C open with O_EXCL and O_CREAT. Linux man page says this conforms to POSIX. If this is not atomic on Windows and files become unreadable, just delete and run again.
        try:
            with open(temp_name, 'x'):
                if not os.path.exists(final_name):
                    can_start = True
                return return_value_func()
        except:
            pass
        time.sleep(1)
        still_blocking = should_block and (datetime.now() - block_start_time).total_seconds() < max_block_time
    
    return return_value_func()


def chunk_finish(final_name):
    temp_name = final_name + '.tmp'
    try:
        os.remove(temp_name)
    except:
        pass


def chunk_start_xarray(base_name, variable_name, should_block=False, max_block_time=float('inf')):
    final_name = xru.get_xarray_filename(base_name, variable_name)
    return chunk_start(final_name, should_block=should_block, max_block_time=max_block_time)

def chunk_finish_xarray(base_name, variable_name):
    final_name = xru.get_xarray_filename(base_name, variable_name)
    temp_name = final_name + '.tmp'
    try:
        os.remove(temp_name)
    except:
        pass


# should_test_chunk_start = False
should_test_chunk_start = True


def debug_raise(given_analysis_filename=None, given_can_start=None):
    if given_analysis_filename is not None and given_can_start == True:
        chunk_finish(given_analysis_filename)
    raise Exception('Debug')
    
def debug_set_trace(given_analysis_filename=None, given_can_start=None):
    if given_analysis_filename is not None and given_can_start == True:
        chunk_finish(given_analysis_filename)
    import pdb
    pdb.set_trace()





if __name__ == '__main__':
    
    # Testing
    
    if should_test_chunk_start:
        # test_chunk_start_silent = False
        test_chunk_start_silent = True
        
        can_start = None
        final_exists = None
        temp_name = None
        test_fh, test_base_name = tempfile.mkstemp(prefix='atomic_test', dir='.')
        os.close(test_fh)
        os.remove(test_base_name)
        
        def test_print_heading(given_str):
            if test_chunk_start_silent: return
            print(); print(given_str)
        
        def test_print_state():
            if test_chunk_start_silent: return
            print(f'can_start {can_start}, final_exists {final_exists}, temp_name {temp_name}')
            print('Listing relevant files')
            print('\n'.join(glob.glob(f'{test_base_name}*')))
        
        test_print_heading('Test chunk_start')
        
        test_final_name = f'{test_base_name}.nc'
        
        test_print_heading('Trying first time')
        can_start, final_exists, temp_name = chunk_start(test_final_name)
        test_print_state()
        assert can_start
        assert not final_exists
        assert temp_name == test_final_name + '.tmp'
        assert os.path.exists(temp_name)
        assert not os.path.exists(test_final_name)
        
        if can_start:
            test_print_heading('Can start, saving')
            xru.reformat_xarray_for_saving(xr.DataArray([np.pi])).to_netcdf(test_final_name)
            chunk_finish(test_final_name);
            test_print_state()
            assert not os.path.exists(temp_name)
            assert os.path.exists(test_final_name)
        
        test_print_heading('Trying second time')
        can_start, final_exists, temp_name = chunk_start(test_final_name)
        test_print_state()
        assert not os.path.exists(temp_name)
        assert os.path.exists(test_final_name)
        
        test_print_heading('Removing test file')
        os.remove(test_final_name)
        test_print_state()
        assert not os.path.exists(test_final_name)
        
        
        test_print_heading('Test chunk_start_xarray')
        
        test_variable_name = 'b'
        test_final_name = xru.get_xarray_filename(test_base_name, test_variable_name)
        
        can_start = None
        final_exists = None
        temp_name = None
        
        test_print_heading('Trying first time')
        # can_start, final_exists, temp_name = chunk_start(test_final_name)
        can_start, final_exists, temp_name = xru.chunk_start_xarray(test_base_name, test_variable_name)
        test_print_state()
        assert can_start
        assert not final_exists
        assert temp_name == test_final_name + '.tmp'
        assert os.path.exists(temp_name)
        assert not os.path.exists(test_final_name)
        
        if can_start:
            test_print_heading('Can start, saving')
            # xru.save_xarray(test_final_name, a=xr.DataArray([np.pi]))
            xru.save_xarray(test_base_name, **{test_variable_name: xr.DataArray([np.pi])})
            # chunk_finish(test_final_name);
            xru.chunk_finish_xarray(test_base_name, test_variable_name);
            test_print_state()
            assert not os.path.exists(temp_name)
            assert os.path.exists(test_final_name)
        
        test_print_heading('Trying second time')
        # can_start, final_exists, temp_name = chunk_start(test_final_name)
        can_start, final_exists, temp_name = xru.chunk_start_xarray(test_base_name, test_variable_name)
        test_print_state()
        assert not os.path.exists(temp_name)
        assert os.path.exists(test_final_name)
                
        test_print_heading('Removing test file')
        os.remove(test_final_name)
        test_print_state()
        assert not os.path.exists(test_final_name)

        # raise NotImplementedError
        # debug_raise()
    
