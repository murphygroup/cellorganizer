/*
 * Copyright (C) 2006 Murphy Lab,Carnegie Mellon University
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published
 * by the Free Software Foundation; either version 2 of the License,
 * or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA.
 * 
 * For additional information visit http://murphylab.web.cmu.edu or
 * send email to murphy@cmu.edu
 */
////////////////////////////////////////////////////////////////////////////////
//  Calculates moments for ALL objects in an image in one pass; replacement for
//  mb_imgmoments called in mb_imgfeatures, where mb_imgmoments calculates
//  moments for one object per call (requires an image containing 1 object only
//
//  [M00, M10, M01, OBJS] = gp_moments_1(IMAGE, LABELLED_IMAGE)
//  where:
//     -M00, M10, M01 are 1xnum_of_objects matrices containing moments for each
//      object
//     -OBJS is a 1xnum_of_objects matrix containing the size of each object
//     -IMAGE is an image of type UINT8 containing N objects
//     -LABELLED_IMAGE is a matrix with size==IMAGE, of type INT32, generated
//      by BWLABEL(IMAGE) after binarization of IMAGE
//
//   08-Aug-01  G. Porreca
//
////////////////////////////////////////////////////////////////////////////////

#include "mex.h"

void mexFunction(int nlhs,               //number of pointers to return args 
		 mxArray *plhs[],        //vector of pointers to return args 
		 int nrhs,               //number of pointers to input  args
		 const mxArray *prhs[]){ //vector of pointers to input  args

  //INPUT ARGUMENTS
  int *image;          // original image
  int *labeled;        // labeled object image

  //OUTPUT ARGUMENTS
  double *a;            // moment00
  double *b;            // moment10 (X)
  double *c;            // moment01 (Y)
  double *d;            // object sizes

  //LOCAL VARIABLES
  int f, index, yrows, xcols, moment_length ;
  


  //perform type checking on input matrices
  if( !mxIsInt32( prhs[0])) 
    mexErrMsgTxt("Input IMAGE must be of class int32");
  image = (int*) mxGetData( prhs[0]);

  if( !mxIsInt32( prhs[1])) 
    mexErrMsgTxt("Input LABELLED IMAGE must be of class int32");
  labeled = (int*) mxGetData( prhs[1]);
  

  //size of input image (and labeled image)
  yrows = mxGetM(prhs[0]);
  xcols = mxGetN(prhs[0]);


  //size of output moment matrices (num_of_objects)
  moment_length = number_of_objects(labeled, xcols, yrows);


  //assign mxArrays to output arguments
  plhs[0] = mxCreateDoubleMatrix(1,moment_length, mxREAL);
  plhs[1] = mxCreateDoubleMatrix(1,moment_length, mxREAL);
  plhs[2] = mxCreateDoubleMatrix(1,moment_length, mxREAL);
  plhs[3] = mxCreateDoubleMatrix(1,moment_length, mxREAL);
  

  //get C pointers to output arguments
  a = mxGetPr(plhs[0]);
  b = mxGetPr(plhs[1]);
  c = mxGetPr(plhs[2]);
  d = mxGetPr(plhs[3]);
  

  //calculate moments and object sizes
  f = calc_Moments(image, labeled, xcols, yrows, moment_length, a, b, c, d);
} //end mexFunction(...)





// determine the number of objects in the labeled image; i.e. the object 
// with the highest number + 1 (object[0:N-1], num_of_objs[1:N])
int number_of_objects(int *label,     // labeled object image 
		      int xcols,      // number of columns in image
		      int yrows){     // number of rows in image

  int i, j, index=0, high;


  for(i = 0; i < xcols; i++){
    for(j = 0; j < yrows; j++){

      if( *(label + index) > high){
	high = *(label + index);
      }

      index++;
    } //end for j
  } //end for i

  return high + 1;
} //end number_of_objects(...)





// calculate the moments for each object in the image in one pass
int calc_Moments(int *img,            // original image
		 int *labl,           // labeled object image
		 int xcols,           // number of columns in img
		 int yrows,           // number of rows in img
		 int num_objs,        // number of objects (1:N) in img
		 double *a,           // moment00 
		 double *b,           // moment10 (X)
		 double *c,           // moment01 (Y)
		 double *d){          // object sizes

  int i, index=0, j, moment_index;

  //return arrays will be incremented below; initialize them
  for(i = 0; i < num_objs; i++){
    *(a + i) = 0;                     
    *(b + i) = 0;                     
    *(c + i) = 0;                     
    *(d + i) = 0;                      
  } //end for i


  //iterate on the number of pixels
  for(i = 0; i < xcols; i++){
    for(j = 0; j < yrows; j++){

      index = (i * yrows) + j;
      moment_index = *(labl + index);

      //sum of fluorescence for the current object
      *(a + moment_index) = *(a + moment_index) + *(img + index);

      //X-weighted sum of fluorescence (i) for the current object
      *(b + moment_index) = *(b + moment_index) + (*(img + index) * (i+1));

      //Y-weighted sum of fluorescence (j) for the current object
      *(c + moment_index) = *(c + moment_index) + (*(img + index) * (j+1));

      //increment the number of pixels (size) for the current object
      *(d + moment_index) = *(d + moment_index) + 1;

    } //end for j
  } //end for i

  return 0;
} //end calcMoments(...) 
						         





 












