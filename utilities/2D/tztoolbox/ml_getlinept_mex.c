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

/*Get all points on a line
 * pts = ml_getlinept_mex(s,t)
 * output: an nx2 array
 * input: 1x2 array for start point and 1x2 array for ending point. s should be
 *     less that t for both x and y
 */

#include "mex.h"
#include "stdio.h"

void mexFunction(int nlhs,               
		 mxArray *plhs[],        
		 int nrhs,               
		 const mxArray *prhs[]){ 
  double *ds,*dt;
  int s[2],t[2];
  int dx,dy,stepx,stepy,maxsize;
  double *pts;
  int npts = 0;
  int fraction;

  ds = (double *)mxGetData(prhs[0]);
  dt = (double *)mxGetData(prhs[1]);

  s[0] = (int)ds[0]; s[1] = (int)ds[1];
  t[0] = (int)dt[0]; t[1] = (int)dt[1];

  dx=t[0]-s[0];
  dy=t[1]-s[1];
  
  maxsize = dy+dx+1;

  plhs[0] = mxCreateDoubleMatrix(2,maxsize,mxREAL);
  pts = mxGetPr(plhs[0]);

  if (dy < 0) {
    dy = -dy;
    stepy = -1;
  }
  else 
    stepy = 1;
    
  if (dx < 0) {
    dx = -dx;  
    stepx = -1;
  }
  else 
    stepx = 1;
    
  dy=dy*2;
  dx=dx*2;
  
  npts++;
  *(pts++) = s[0];
  *(pts++) = s[1];
    
  if (dx > dy) {
    fraction = dy - dx/2;
    while(s[0]!=t[0]) {
      if(fraction >= 0) {
	s[1]=s[1]+stepy;
	fraction=fraction-dx;
      }
      s[0]=s[0]+stepx;
      fraction=fraction+dy;
      npts++;
      *(pts++) = s[0];
      *(pts++) = s[1];
    }
  }
  else {
    fraction= dx - dy/2;
    while (s[1] != t[1]) {
      if (fraction >= 0) {
	s[0]=s[0]+stepx;
	fraction=fraction-dy;
      }
      s[1] = s[1]+stepy;
      fraction=fraction+dx;
      npts++;
      *(pts++) = s[0];
      *(pts++) = s[1];
    }
  }
  
//   
//   if(nlhs==2) {
//     plhs[1] = mxCreateDoubleMatrix(1,1,mxREAL);
//     double* pnpts = mxGetPr(plhs[1]);
//     *(pnpts) = npts;
//   }
    
}
