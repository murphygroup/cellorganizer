/*/////////////////////////////////////////////////////////////////////////
//
//
//                            ml_3Dcoocmat.c 
//
//
//                            Xiang Chen
//                            modified from mb_texture.c
//                            Aug 19, 2002
//
//  Revisions:
//  modified October_2008: ml_3Dcoocmat.c return only the cooccurrence matrix
//
/////////////////////////////////////////////////////////////////////////*/

#include "mex.h"
#include "math.h"
#include "matrix.h"
#include "Include/ppgm.h"
#include "Include/3DCVIPtexture.h"
#include "sys/types.h"
#include "stdlib.h"
#include "memory.h"

#define Y 0
#define X 1
#define Z 2

#define row 0
#define col 1

#define ind(X, Y, Z)             (Y) + (X) * ny + (Z) * ny * nx


/*
#define SIGN(X,Y) ((Y)<0 ? -fabs(X) : fabs(X))
#define DOT fprintf(stderr,".")
#define SWAP(a,b) {Y=(a);(a)=(b);(b)=Y;}
*/

#define idx(X, Y, Z) (Y) + (X) * ny + (Z) * ny * nx
#define TRUE  1
#define FALSE 0

float **pgm_matrix ();



void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
  int         distance;                 /*parameter for texture calculations*/
  u_int8_t*   p_img;                    /*The image from Matlab*/
  int         ny;                       /*Image y*/
  int         nx;                       /*Image x*/
  int         nz;                       /*Image z*/
  int         NDims;
  long        Nvoxels;
  const int*  dims;
  int         Dims[3];
  int         i, j, k,x,y,z;
  long        offset ;                  
  int         outputsize[2] ;           /*Dimensions of TEXTURE struct*/
  float*      output ;
  float**     P_matrix;					/* cooccurrence matrix */
  
  int*	      grayLevels;				/* vector of gray levels */
  int*	      tmpGrayLevels;
  int 		  nbGL;                     /* number of gray levels */
  int         NDimDir;
  int         dimDir;

  int**       directions;				/* Direction/Distance to take into account, each row has 3 values to define x,y,z distances */	
  int*        direct;
  
  int tonec[PGM_MAXMAXVAL+1];
  int R[13];
  int directionToCount[13];
  int argn, bps, padright, m, d;
  int itone, jtone, g_val,tones;
  int row_,col_;
  
  gray nmaxval, maxval;

  d = 1;

/*
nrhs: number of input 
nlhs: number of output
prhs[]: array of input
plhs[]: array of output
*/



  if (nrhs != 2) {
    mexErrMsgTxt("ml_3Dcoocmat requires 2 input arguments.\n") ;
  } else if (nlhs != 2) {
    mexErrMsgTxt("ml_3Dcoocmat returns 2 output arguments.\n") ;
  }

  if (!mxIsNumeric(prhs[0])) {
    mexErrMsgTxt("ml_3Dcoocmat requires nbGL sinbGLe numeric input as nbGL first argument.\n") ;
  }
  if (!mxIsUint8(prhs[0])) {
    mexErrMsgTxt("ml_3Dcoocmat requires nbGL sinbGLe input of type unsigned 8-bit integer.\n") ;
  }
  if (!mxIsNumeric(prhs[1])) {
    mexErrMsgTxt("ml_3Dcoocmat requires nbGL numeric input as nbGL second argument.\n") ;
  }
  if (!mxIsInt32(prhs[1])) {
    mexErrMsgTxt("ml_3Dcoocmat requires as second input nbGL matrix of type int32 .\n") ;
  }

  NDims = mxGetNumberOfDimensions(prhs[0]);
  if (NDims != 3) {
    mexErrMsgTxt("ml_3Dcoocmat requires nbGL 3D image as the 1st input.\n");
  }
  
  NDimDir = mxGetNumberOfDimensions(prhs[1]);
  if (NDimDir > 2) {
    mexErrMsgTxt("ml_3Dcoocmat requires nbGL 2D matrix as the 2nd input.\n");
  }
  dimDir  = mxGetN(prhs[1]);
  if (dimDir != 3) {
    mexErrMsgTxt("ml_3Dcoocmat requires nbGL 2D maxtrix as the 2nd input with 3 columns for x,y,z directions to connsider.\n");
  }

  
  dims = mxGetDimensions(prhs[0]);
  Dims[0] = dims[0];
  Dims[1] = dims[1];
  Dims[2] = dims[2];
  ny = Dims[0];
  nx = Dims[1];
  nz = Dims[2];

  /*mexPrintf("ny=Dims[0]=%i nx=Dims[1]=%i nz=Dims[2]=%i\n",ny,nx,nz);*/

  if(!(nx > 1) || !(ny > 1) || !(nz > 1)) {
    mexErrMsgTxt("ml_3Dcoocmat requires an input 3D image, not nbGL scalar.\n") ;
  }

  Nvoxels = nx * ny * nz;
  /*mexPrintf("Nvoxels=%i\n",Nvoxels);*/

  /* p_img: monodimentional array reprensenting the image along x, y, and z*/
  p_img = (u_int8_t*)mxGetData(prhs[0]) ;

  distance = 1 ;

  /* Get the different direction/distance*/
  dimDir = mxGetM(prhs[1]);
  direct = (int*)mxGetData(prhs[1]) ;


  directions = mxCalloc(dimDir, sizeof(int*)) ;
  if(directions) {
  	for(i=0; i<dimDir ; i++) {
      directions[i] = mxCalloc(3, sizeof(int)) ;
	  if(!directions[i]) mexErrMsgTxt("ml_3Dcoocmat : error allocating directions[i]") ;
    }
  } else mexErrMsgTxt("ml_3Dcoocmat : error allocating directions");
 
  /* directions: fill with the values from the monodimentional array direct
  */
  for(i=0; i<dimDir; i++){ 
    for (j=0;j<3;j++) { 
	   offset = j*dimDir+i;
	   directions[i][j] = (int)direct[offset];
	}
	/*mexPrintf("%i, direction[%i %i %i]\n",i,directions[i][0],directions[i][1],directions[i][2]);*/
  }
    
  
  /**************************************************************************/
  /* CALCULATE THE COOCCURRENCE MATRIX GIVEN THE DIRECTION                  */
  /* all the directions are merged in nbGL*nbGL unique Cooccurrence matrix  */
  /**************************************************************************/
 
 /*int* Calculate_Cooccurrence_Matrix( register gray *p_img, int nx, int ny, int nz, int** directions, int dimDir, float** P_matrix, int *nbGraylevels,float* output, int *ab) 
 grayLevels = Calculate_Cooccurrence_Matrix( p_img,nx,ny,nz,directions,dimDir,P_matrix,&nbGL,&output,add);  
  
  */

   /* Determine the number of different gray scales (not maxval) */
  for (row_ = PGM_MAXMAXVAL; row_ >= 0; --row_)
    tonec[row_] = -1;
  for (k = 0; k < nz; ++k)
    for (i = nx - 1; i >= 0; --i)
      for (j = 0; j < ny; ++j)
	{
	  /*   if (p_img[row_][col_])   If gray value equal 0 don't include */		
	  tonec[p_img[idx(i, j, k)]] = p_img[idx(i, j, k)];
      }	
  
  for (row_ = PGM_MAXMAXVAL, tones = 0; row_ >= 0; --row_)
    if (tonec[row_] != -1)
      tones++;
  /*fprintf (stderr, "(Image has %d graylevels.)\n", tones); */

  nbGL = tones;  
  grayLevels = (int*) malloc(tones*sizeof(int));
 
  /* Collapse array, taking out all zero values */
  i=0;
  for (row_ = 0, itone = 0; row_ <= PGM_MAXMAXVAL; row_++)
    if (tonec[row_] != -1) {
      tonec[row_] = itone++; /* convertion table*/
      *(grayLevels+i) = row_;
      i++;
    }
  /* Now array contains only the gray levels present (in ascending order) */
  

  /* Allocate memory for gray-tone spatial dependence matrix */
  P_matrix = pgm_matrix (0, tones-1, 0, tones-1);
  for (row_ = 0; row_ < tones; ++row_)
    for (col_ = 0; col_ < tones; ++col_)
	  P_matrix[row_][col_] = 0;
	
  for (i = 0; i < 13; i++) {
    R[i] = 0;
    directionToCount[i]= FALSE;
  }

 /* Find gray-tone spatial dependence matrix */
  /* XC: From (x, y, z):
     R[0]:  (x + 1, y, z);
	 R[1]:  (x, y + 1, z);
	 R[2]:  (x + 1, y + 1, z);
	 R[3]:  (x + 1, y - 1, z);
	 R[4]:  (x, y, z + 1);
	 R[5]:  (x + 1, y, z + 1);
	 R[6]:  (x, y + 1, z + 1);
	 R[7]:  (x + 1, y + 1, z + 1);
	 R[8]:  (x + 1, y - 1, z + 1);
	 R[9]:  (x + 1, y, z - 1);
	 R[10]: (x, y + 1, z - 1);
	 R[11]: (x + 1, y + 1, z - 1);
	 R[12]: (x + 1, y - 1, z - 1).
  */
  
  for(m=0; m<dimDir; m++){
  	i = directions[m][0];
  	j = directions[m][1];
  	k = directions[m][2];
  	if (i==0){
  		if (j==0){
  			if (k==1)
  				directionToCount[4]=TRUE;
  		}
  		if (j==1){
  			if (k==0)
  				directionToCount[1]=TRUE;
  			if (k==1)
  				directionToCount[6]=TRUE;
  			if (k==-1)
  				directionToCount[10]=TRUE;	
  		}
  	}else{ /* if (i==1) */
  		if (j==0){
  			if (k==0)
  				directionToCount[0]=TRUE;
  			if (k==1)
  				directionToCount[5]=TRUE;
  			if (k==-1)
  				directionToCount[9]=TRUE;
  		}
  		if (j==1){
  			if (k==0)
  				directionToCount[2]=TRUE;
  			if (k==1)
  				directionToCount[7]=TRUE;	
  			if (k==-1)
  				directionToCount[11]=TRUE;		
  		}
  		if (j==-1){
  			if (k==0)
  				directionToCount[3]=TRUE;
  			if (k==1)
  				directionToCount[8]=TRUE;
  			if (k==-1)
  				directionToCount[12]=TRUE;	
  		}
  	}		
  }

  
 /* fprintf (stderr, "(Computing spatial dependence matrix..."); */
  for (k = 0; k < nz; k++)
    for (j = 0; j < ny; j++)
      for (i = 0; i < nx; i++) {
		/* if (p_img[idx(i, j, k)]) { */
	  	x = tonec[p_img[idx(i, j, k)]];
	  if (directionToCount[0] && i + d < nx) {
	    	y = tonec[p_img[idx(i + d, j, k)]];
	    	P_matrix[x][y]++;
	    	P_matrix[y][x]++;
	    	R[0]+=2;
	  } 
	  if (directionToCount[1] && j + d < ny) {
	    y = tonec[p_img[idx(i, j + d, k)]];
	    P_matrix[x][y]++;
	    P_matrix[y][x]++;
	    R[1]+=2;
	  }
	  if (directionToCount[2] && i + d < nx && j + d < ny) {
	    y = tonec[p_img[idx(i + d, j + d, k)]];
	    P_matrix[x][y]++;
	    P_matrix[y][x]++;
	    R[2]+=2;
	  }
	  if (directionToCount[3] && i + d < nx && j - d >= 0) {
	    y = tonec[p_img[idx(i + d, j - d, k)]];
	    P_matrix[x][y]++;
	    P_matrix[y][x]++;
	    R[3]+=2;
	  } 
	  if (directionToCount[4] && k + d < nz) {
	    y = tonec[p_img[idx(i, j, k + d)]];
	    P_matrix[x][y]++;
	    P_matrix[y][x]++;
	    R[4]+=2;
	  } 
	  if (directionToCount[5] && i + d < nx && k + d < nz) {
	    y = tonec[p_img[idx(i + d, j, k + d)]];
	    P_matrix[x][y]++;
	    P_matrix[y][x]++;
	    R[5]+=2; 
	  }
	  if (directionToCount[6] && j + d < ny && k + d < nz) {
	    y = tonec[p_img[idx(i, j + d, k + d)]];
	    P_matrix[x][y]++;
	    P_matrix[y][x]++;
	    R[6]+=2;
	  } 
	  if (directionToCount[7] && i + d < nx && j + d < ny && k + d < nz) {
	    y = tonec[p_img[idx(i + d, j + d, k + d)]];
	    P_matrix[x][y]++;
	    P_matrix[y][x]++;
	    R[7]+=2;
	  }
	  if (directionToCount[8] && i + d < nx && j - d >= 0 && k + d < nz) {
	    y = tonec[p_img[idx(i + d, j - d, k + d)]];
	    P_matrix[x][y]++;
	    P_matrix[y][x]++;
	    R[8]+=2;
	  }
	  if (directionToCount[9] && i + d < nx && k - d >= 0) {
	    y = tonec[p_img[idx(i + d, j, k - d)]];
	    P_matrix[x][y]++;
	    P_matrix[y][x]++;
	    R[9]+=2;
	  }
	  if (directionToCount[10] && j + d < ny && k - d >= 0) {
	    y = tonec[p_img[idx(i, j + d, k - d)]];
	    P_matrix[x][y]++;
	    P_matrix[y][x]++;
	    R[10]+=2;
	  }
	  if (directionToCount[11] && i + d < nx && j + d < ny && k - d >= 0) {
	    y = tonec[p_img[idx(i + d, j + d, k - d)]];
	    P_matrix[x][y]++;
	    P_matrix[y][x]++;
	    R[11]+=2;
	  }
	  if (directionToCount[12] && i + d < nx && j - d >= 0 && k - d >= 0) {
	    y = tonec[p_img[idx(i + d, j - d, k - d)]];
	    P_matrix[x][y]++;
	    P_matrix[y][x]++;
	    R[12]+=2;
	  }
/*	} */
      }
  /* Gray-tone spatial dependence matrices are complete */
  
  
  /* Normalize gray-tone spatial dependence matrix */
/*  for (k = 1; k < 13; k++) {
  	R[0]+=R[k];
  }    
  
  for (itone = 0; itone < tones; ++itone)
    for (jtone = 0; jtone < tones; ++jtone)
    	P_matrix[itone][jtone] /= R[0];
*/

 
 /*****************************************/


  outputsize[row] = nbGL;
  outputsize[col] = nbGL;
  
  plhs[0] = mxCreateNumericArray(2, outputsize, mxSINGLE_CLASS, mxREAL) ;
  if (!plhs[0]) mexErrMsgTxt("ml_3Dcoocmat: error allocating return variable.") ;
  output = (float*)mxGetData(plhs[0]) ; 
  
  for(j=0; j<nbGL; j++)
	for (i=0; i<nbGL; i++){
		offset = j*nbGL+i;
		output[offset] = P_matrix[i][j];
		/*mexPrintf("output[%i]=%f\n",offset,output[offset]);*/
	}

 
  outputsize[row]=nbGL;
  outputsize[col]=1;
  plhs[1] = mxCreateNumericArray(2, outputsize, mxINT32_CLASS, mxREAL) ;
  tmpGrayLevels = (int*)mxGetData(plhs[1]) ; 
  for (i=0; i<nbGL; i++){
		tmpGrayLevels[i] = grayLevels[i];
		/*mexPrintf("tmpGrayLevels[%i]=%i\n",i,tmpGrayLevels[i]);*/
	}

  /*
    Memory clean-up.
  */
  for(i=0; i<dimDir ; i++) {
    mxFree(directions[i]);
  }
  mxFree(directions);
  
  /*
  for (i=0; i<nbGL; i++)
    mxFree(P_matrix[i]);
  mxFree(P_matrix);
  mexPrintf("Hello5\n");
 
  mxFree(grayLevels);
  mexPrintf("Hello6\n");
  */
}



float **pgm_matrix (nrl, nrh, ncl, nch)
  int nrl, nrh, ncl, nch;

/* Allocates a float matrix with range [nrl..nrh][ncl..nch] */
{
  int i;
  float **m;

  /* allocate pointers to rows */
  m = (float **) malloc ((unsigned) (nrh - nrl + 1) * sizeof (float *));
  if (!m)
    fprintf (stderr, "memory allocation failure (pgm_matrix 1) "), exit (1);
  m -= ncl;

  /* allocate rows and set pointers to them */
  for (i = nrl; i <= nrh; i++)
  {
    m[i] = (float *) malloc ((unsigned) (nch - ncl + 1) * sizeof (float));
    if (!m[i])
      fprintf (stderr, "memory allocation failure (pgm_matrix 2) "), exit (2);
    m[i] -= ncl;
  }
  /* return pointer to array of pointers to rows */
  return m;
}




