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
# include "mex.h"

void ml_hitmiss(double *changed_image, double *oper, double *struct_elem, double *image_exp, int r, int c)
{

int hits = 0;
int n, m, z, y, x, w;
*oper = 0;

   for (n = 2; n < c+2 ; n++)
   {
      for(m = 2; m < r+2; m++)
      {          
          x = n - 1;      
          for(w = 0; w < 3; w++)
          {      
             y = m - 1;    
             for(z = 0; z < 3; z++)
             {  

                if(struct_elem[w+z*3] == image_exp[x+y*(c+4)]) 
                    hits = hits + 1;             
                else if(struct_elem[w+z*3] == 2)     
                    hits = hits + 1;                                  

               y = y + 1;                
             }
             
             x = x + 1;
          }

          if (hits == 9)
          {
             changed_image[(n-2)+(m-2)*c] = 1;
             *oper = *oper + 1;

          } else {

             changed_image[(n-2)+(m-2)*c] = 0;
	  }

          hits = 0;
      }
   }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

	mxArray *struct_elem_Array, *image_exp_Array;
	double *struct_elem, *image_exp;
	double *changed_image;
	double *oper;
	int r, c;
	int x,y;

	
	struct_elem_Array = prhs[0];
	struct_elem = mxGetPr(struct_elem_Array);

	image_exp_Array = prhs[1];
	r = mxGetN(image_exp_Array)-4;
	c = mxGetM(image_exp_Array)-4;
	image_exp = mxGetPr(image_exp_Array);

	plhs[0] = mxCreateDoubleMatrix(c, r, mxREAL);
	changed_image = mxGetPr(plhs[0]);

	plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
	oper = mxGetPr(plhs[1]);

        ml_hitmiss(changed_image, oper, struct_elem, image_exp, r, c);

}
