#include "./armadillo-6.700.6/mex_interface/armaMex.hpp"
#include "mex.h" 
#include "iostream"
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs != 1)
        mexErrMsgTxt("Incorrect number of input arguments");
    
    if ((mxGetClassID(prhs[0]) != mxDOUBLE_CLASS))
        mexErrMsgTxt("Input must be of type double.");
    
    if ((mxIsComplex(prhs[0]))) 
        mexErrMsgTxt("Input must be real.");
        
    Mat<double> f = armaGetData<double>(prhs[0]);
    //Mat<double> f = f0;
    //cout << f << endl;
   
    int l = f.n_rows;
    int m = f.n_cols;
    Mat<double> g = zeros(l, m);

    g.rows(1, l - 2) = (diff(f.rows(1, l - 1), 1, 0) - diff(f.rows(0, l - 2), 1, 0)) / 2.0;
    g.row(0) = g.row(1) * 2 - g.row(2);
    g.row(l - 1) = g.row(l - 2) * 2 - g.row(l - 3);
    mat v = g;
    //g = g.zeros();
    g.cols(1, m - 2) = (diff(f.cols(1, m - 1), 1, 1) - diff(f.cols(0, m - 2), 1, 1)) / 2.0;
    g.col(0) = g.col(1) * 2 - g.col(2);
    g.col(m - 1) = g.col(m - 2) * 2 - g.col(m - 3);

    v = (v + g) / 2.0;
    plhs[0] = armaCreateMxMatrix(l, m);
    armaSetData(plhs[0], v);
    return;                   
}			
				
				

			
				
	
