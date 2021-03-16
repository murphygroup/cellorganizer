/*
 * Utility functions for IPT MEX-files.
 *
 * Copyright 1993-2002 The MathWorks, Inc.
 */

/* $Revision: 1.10 $ */

#include <math.h>
#include "mex.h"

static char rcsid[] = "$Id: iptutil.c,v 1.10 2002/03/15 15:28:27 eddins Exp $";

/*
 * check_nargin --- interface to IPT checknargin function for 
 * checking for the proper number of input arguments.
 */
void check_nargin(double low, 
                  double high, 
                  int numInputs, 
                  const char *function_name)
{
    mxArray *prhs[4];
    mxArray *plhs[1];
    int nlhs = 0;
    int nrhs = 4;

    plhs[0] = NULL;

    prhs[0] = mxCreateScalarDouble((double) low);
    prhs[1] = mxCreateScalarDouble((double) high);
    prhs[2] = mxCreateScalarDouble((double) numInputs);
    prhs[3] = mxCreateString(function_name);

    mexCallMATLAB(nlhs, plhs, nrhs, prhs, "checknargin");
    
    mxDestroyArray(prhs[0]);
    mxDestroyArray(prhs[1]);
    mxDestroyArray(prhs[2]);
    mxDestroyArray(prhs[3]);

    if (plhs[0] != NULL)
    {
        mxDestroyArray(plhs[0]);
    }
}

/*
 * check_input --- interface to IPT checkarray function for
 * checking validity of input arguments.
 */
void check_input(const mxArray *A,
                 const char    *classes,
                 const char    *attributes,
                 const char    *function_name,
                 const char    *variable_name,
                 int           argument_position)
{
    mxArray *prhs[6];
    mxArray *plhs[1];
    int nlhs = 0;
    int nrhs = 6;
        
    prhs[0] = (mxArray *) A;
    prhs[1] = mxCreateString(classes);
    prhs[2] = mxCreateString(attributes);
    prhs[3] = mxCreateString(function_name);
    prhs[4] = mxCreateString(variable_name);
    prhs[5] = mxCreateScalarDouble((double) argument_position);

    plhs[0] = NULL;
    
    mexCallMATLAB(nlhs, plhs, nrhs, prhs, "checkinput");

    mxDestroyArray(prhs[1]);
    mxDestroyArray(prhs[2]);
    mxDestroyArray(prhs[3]);
    mxDestroyArray(prhs[4]);
    mxDestroyArray(prhs[5]);

    if (plhs[0] != NULL)
    {
        mxDestroyArray(plhs[0]);
    }
}

mxArray *call_one_input_one_output_function(const mxArray *A, 
                                            const char *function_name)
{
    int nrhs = 1;
    int nlhs = 1;
    mxArray *prhs[1];
    mxArray *plhs[1];

    prhs[0] = (mxArray *) A;
    
    mexCallMATLAB(nlhs, plhs, nrhs, prhs, function_name);
    
    return plhs[0];
}
    
mxArray *convert_to_logical(const mxArray *A)
{
    int nrhs = 2;
    int nlhs = 1;
    mxArray *prhs[2];
    mxArray *plhs[1];

    prhs[0] = (mxArray *) A;
    prhs[1] = mxCreateScalarDouble(0.0);
    
    mexCallMATLAB(nlhs, plhs, nrhs, prhs, "ne");

    mxDestroyArray(prhs[1]);
    
    return plhs[0];
}
    
/*
 * get_array_element
 *
 * Given a void pointer to the beginning of an array's data, the class ID
 * of the array, and a linear element offset, return the corresponding 
 * element value as a double.
 * 
 */
double get_array_element(void *ptr, mxClassID array_class, int offset)
{
    double result;

    switch (array_class)
    {
    case mxUINT8_CLASS:
        result = *((uint8_T *) ptr + offset);
        break;

    case mxLOGICAL_CLASS:
        result = *((uint8_T *) ptr + offset);
        break;
        
    case mxUINT16_CLASS:
        result = *((uint16_T *) ptr + offset);
        break;
        
    case mxUINT32_CLASS:
        result = *((uint32_T *) ptr + offset);
        break;
        
    case mxINT8_CLASS:
        result = *((int8_T *) ptr + offset);
        break;
        
    case mxINT16_CLASS:
        result = *((int16_T *) ptr + offset);
        break;
        
    case mxINT32_CLASS:
        result = *((int32_T *) ptr + offset);
        break;
        
    case mxSINGLE_CLASS:
        result = *((float *) ptr + offset);
        break;
        
    case mxDOUBLE_CLASS:
        result = *((double *) ptr + offset);
        break;
        
    default:
        mexErrMsgTxt("Unsupported input class.");
    }
    
    return result;
}


/*
 * Macros to safely convert floating point values to integer ranges.
 * This macro is necessary because not all C compilers do the same
 * thing for NaN's or out-of-range conversions.
 */
#define CONVERT_SIGNED(p,val,type,min,max) \
  if (mxIsNaN(val)) val = 0.0;    \
  if (val > max)    val = max;    \
  if (val < min)    val = min;    \
  if (val < 0.0)    *(((type *) p) + offset) = (type)(val - 0.5); \
  else              *(((type *) p) + offset) = (type)(val + 0.5)
 
#define CONVERT_UNSIGNED(p,val,type,min,max) \
  if (mxIsNaN(val)) val = 0.0;    \
  if (val > max)    val = max;    \
  if (val < min)    val = min;    \
  *(((type *) p) + offset) = (type)(val + 0.5)

/*
 * put_array_element
 *
 * Given a void pointer to the beginning of an array's data, a double-precision
 * value, the class of the array, and a linear element offset, set the
 * corresponding array element to the given value.  If assigning to
 * an integer array, convert NaNs to 0, clip out-of-range values to the
 * target range, and round.
 */
void put_array_element(void *pr_a, double val, mxClassID class_a,
                       int offset)
{
    switch (class_a)
    {
    case mxDOUBLE_CLASS:
        *((double *) pr_a + offset) = val;
        break;
        
    case mxSINGLE_CLASS:
        *((float *) pr_a + offset) = (float) val;
        break;
        
    case mxUINT8_CLASS:
        CONVERT_UNSIGNED(pr_a, val, uint8_T, MIN_uint8_T, MAX_uint8_T);
        break;
        
    case mxLOGICAL_CLASS:
        CONVERT_UNSIGNED(pr_a, val, uint8_T, MIN_uint8_T, MAX_uint8_T);
        break;
        
    case mxUINT16_CLASS:
        CONVERT_UNSIGNED(pr_a, val, uint16_T, MIN_uint16_T, MAX_uint16_T);
        break;
        
    case mxUINT32_CLASS:
        CONVERT_UNSIGNED(pr_a, val, uint32_T, MIN_uint32_T, MAX_uint32_T);
        break;
        
    case mxINT8_CLASS:
        CONVERT_SIGNED(pr_a, val, int8_T, MIN_int8_T, MAX_int8_T);
        break;
        
    case mxINT16_CLASS:
        CONVERT_SIGNED(pr_a, val, int16_T, MIN_int16_T, MAX_int16_T);
        break;
        
    case mxINT32_CLASS:
        CONVERT_SIGNED(pr_a, val, int32_T, MIN_int32_T, MAX_int32_T);
        break;
        
    default:
        mexErrMsgTxt("Unsupported input class.");
    }
}

/*
 * Given a linear offset p into an M-by-N array, row and column offsets 
 * r and c, and matrix dimensions M and N, return false if the 
 * corresponding neighbor is outside the bounds of the array; 
 * otherwise return true.
 */
bool is_inside(int p, int r_offset, int c_offset, int M, int N)
{
    int row = p % M;
    int col = p / M;
    int new_col = col + c_offset;
    int new_row = row + r_offset;
    
    return ((new_col >= 0) && (new_col < N) &&
            (new_row >= 0) && (new_row < M));
}

/*
 * Given input style (either 4 or 8) and row dimension M, 
 * initialize neighbor arrays nhood_r (row offsets), nhood_c 
 * (column offsets), and nhood (linear offsets).  Also initialize 
 * num_neighbors.
 *
 * nhood_r, nhood_c, and nhood need to have room for at least
 * eight ints.
 */
extern void init_neighbors(int style, int M, int nhood_r[], int nhood_c[], 
                           int nhood[], int *num_neighbors)
{
    int k;
    
    switch (style)
    {
    case 8:
        /* N */
        nhood_r[0] = -1;
        nhood_c[0] = 0;

        /* NE */
        nhood_r[1] = -1;
        nhood_c[1] = 1;
        
        /* E */
        nhood_r[2] = 0;
        nhood_c[2] = 1;
        
        /* SE */
        nhood_r[3] = 1;
        nhood_c[3] = 1;
        
        /* S */
        nhood_r[4] = 1;
        nhood_c[4] = 0;
        
        /* SW */
        nhood_r[5] = 1;
        nhood_c[5] = -1;
        
        /* W */
        nhood_r[6] = 0;
        nhood_c[6] = -1;
        
        /* NW */
        nhood_r[7] = -1;
        nhood_c[7] = -1;

        *num_neighbors = 8;
        break;
        
    case 4:
        /* N */
        nhood_r[0] = -1;
        nhood_c[0] = 0;
        
        /* E */
        nhood_r[1] = 0;
        nhood_c[1] = 1;
        
        /* S */
        nhood_r[2] = 1;
        nhood_c[2] = 0;
        
        /* W */
        nhood_r[3] = 0;
        nhood_c[3] = -1;
        
        *num_neighbors = 4;
        break;
        
    default:
        mexErrMsgTxt("Internal problem: invalid style");
    }

    for (k = 0; k < *num_neighbors; k++)
    {
        nhood[k] = M*nhood_c[k] + nhood_r[k];
    }
}


