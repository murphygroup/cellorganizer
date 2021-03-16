/*/////////////////////////////////////////////////////////////////////////
//
//
//                            ml_Znl.cpp 
//
//
//                           Michael Boland
//                            09 Dec 1998
//
//  Revisions: EJSR - ~<complex.h> -> <complex>
//                    +using namespace std;
//                    Name from mb_Znl.cpp
//
/////////////////////////////////////////////////////////////////////////*/


#include "mex.h"
#include "matrix.h"
#include <complex>
#include <math.h>
#include <sys/types.h>

#define row 0
#define col 1

using namespace std;

//
// Calculates n! (uses double arithmetic to avoid overflow)
//
double factorial(double n)
{
	if(n < 0)
		return(0.0) ;
	if(n == 0.0)
		return(1.0) ;
	else
		return(n * factorial(n-1.0)) ;
}	


void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{

  int n ;                        /* degree (n) of the Zernike moment */
  int l ;                        /* angular dependence of the Zernike moment */
  double* X ;                    /* list of X coordinates of pixels */
  double* Y ;                    /* list of Y coordinates of pixels */
  double* P ;                    /* list of values of pixels */
  int outputsize[2] = {1,1} ;    /* Dimensions of return variable */
  register double x, y, p ;      /* individual values of X, Y, P */
  int i ;
  complex<double> sum = 0.0 ;    /* Accumulator for complex moments */
  complex<double> Vnl = 0.0 ;    /* Inner sum in Zernike calculations */
  double* preal ;                /* Real part of return value */
  double* pimag ;                /* Imag part of return value */
  int m ;



  if (nrhs != 5) {
    mexErrMsgTxt("ml_Znl(N, L, X, Y, P), Zernike moment generating function. "
                 "The moment of degree n and angular dependence l for the "
                 "pixels defined by coordinate vectors X and Y and intensity "
                 "vector P.  X, Y, and P must have the same length.") ;

  } else if (nlhs != 1) {
    mexErrMsgTxt("ml_Znl returns a single output.\n") ;
  }

  if ( !mxIsNumeric(prhs[0]) || (mxGetM(prhs[0]) != 1) || 
      (mxGetN(prhs[0]) != 1) ) {
    mexErrMsgTxt("The first argument (n) should be a scalar\n") ;
  }

  if ( !mxIsNumeric(prhs[1]) || (mxGetM(prhs[1]) != 1) || 
      (mxGetN(prhs[1]) != 1) ) {
    mexErrMsgTxt("The second argument (l) should be a scalar\n") ;
  }

  if ( !mxIsNumeric(prhs[2]) || (mxIsComplex(prhs[2])) ) {
    mexErrMsgTxt("The third argument (X) should be numeric and not complex.") ;
  }

  if ( !mxIsNumeric(prhs[3]) || (mxIsComplex(prhs[3])) ) {
    mexErrMsgTxt("The third argument (Y) should be numeric and not complex.") ;
  }

  if ( !mxIsNumeric(prhs[4]) || (mxIsComplex(prhs[4])) ) {
    mexErrMsgTxt("The third argument (P) should be numeric and not complex.") ;
  }

  if (mxGetM(prhs[2])!=mxGetM(prhs[3]) || (mxGetM(prhs[3])!=mxGetM(prhs[4]))){
    mexErrMsgTxt("X, Y, and P must have the same number of rows.") ;
  }

  if (mxGetM(prhs[2]) < mxGetM(prhs[2])) {
    mexErrMsgTxt("X, Y, and P should be column vectors.") ;
  }

  n = (int)mxGetScalar(prhs[0]) ;
  l = (int)mxGetScalar(prhs[1]) ;

  X = mxGetPr(prhs[2]) ;
  Y = mxGetPr(prhs[3]) ;
  P = mxGetPr(prhs[4]) ;

  for(sum = 0.0, i = 0 ; i < mxGetM(prhs[2]) ; i++) {
    x = X[i] ;
    y = Y[i] ;
    p = P[i] ;
    
    for(Vnl = 0.0, m = 0; m <= (n-l)/2; m++) {
      Vnl += (pow((double)-1.0,(double)m)) * ( factorial(n-m) ) / 
	( factorial(m) * (factorial((n - 2.0*m + l) / 2.0)) *
	  (factorial((n - 2.0*m - l) / 2.0)) ) *
	( pow( sqrt(x*x + y*y), (double)(n - 2*m)) ) *
	polar(1.0, l*atan2(y,x)) ;
      //
      // NOTE: This function did not work with the following:
      //  ...pow((x*x + y*y), (double)(n/2 -m))...
      //  perhaps pow does not work properly with a non-integer
      //  second argument.
      // 'not work' means that the output did not match the 'old'
      //  Zernike calculation routines.
      //
    }

    sum += p * conj(Vnl) ;
  }

  sum *= (n+1)/PI ;
  

  /* Assign the returned value */

  plhs[0] = mxCreateNumericArray(2, outputsize, mxDOUBLE_CLASS, mxCOMPLEX) ;
  preal = mxGetPr(plhs[0]) ;
  pimag = mxGetPi(plhs[0]) ;

  preal[0] = real(sum) ;
  pimag[0] = imag(sum) ;

}
