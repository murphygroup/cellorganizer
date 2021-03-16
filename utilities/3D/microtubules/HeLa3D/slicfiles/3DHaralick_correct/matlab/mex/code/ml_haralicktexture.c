/*/////////////////////////////////////////////////////////////////////////
//
//
//                            ml_haralicktexture.c 
//
//
//                            Xiang Chen
//                            modified from mb_texture.c
//                            Aug 19, 2002
//
//  Revisions:
//      May 14, 2003: change name to ml_3Dtexture.c
//      Nov 18, 2008: change name to ml_haralicktexture.c : calculate the 14 features
//                    from the cooccurrence given in parameter 
//
/////////////////////////////////////////////////////////////////////////*/


#include "mex.h"
#include "matrix.h"
#include "Include/ppgm.h"
#include "Include/3DCVIPtexture.h"
#include <sys/types.h>
#include <stdlib.h>
#include <memory.h>

#define y 0
#define x 1
#define z 2

#define row 0
#define col 1

#define ind(x, y, z)   (y) + (x) * ny + (z) * ny * nx


extern TEXTURE * Extract_Texture_Features();

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{

  TEXTURE_FEATURE_MAP* features_used ;  /*Indicate which features to calc.*/
  TEXTURE*    features ;                /*Returned struct of features*/
  char**	  feature_names;
  int         NDims;
  float* 	  coomat;   
  int* 	      grayLevels;
  int*		  nbGrayLevels;
  int*		  tmpGrayLevels;
  
  const int*  dims;
  int         i, j, k ;

  long        offset ;                  
  int         outputsize[1] ;           /*Dimensions of TEXTURE struct*/
  int         outputindex[2] ;
  float*      output ;                  /*Features to return*/


  if (nrhs != 3) {
    mexErrMsgTxt("ml_haralicktexture requires 3 input arguments.\n") ;
  } else if (nlhs != 2) {
    mexErrMsgTxt("ml_haralicktexture returns two outputs.\n") ;
  }

  if (!mxIsNumeric(prhs[0])) {
    mexErrMsgTxt("ml_haralicktexture requires a single numeric input.\n") ;
  }

  if (!mxIsSingle(prhs[0])) {
    mexErrMsgTxt("ml_haralicktexture requires a single input of type SINGLE (single precision floating point).\n") ;
  }
  if (!mxIsNumeric(prhs[1])) {
    mexErrMsgTxt("ml_haralicktexture requires a vector of gray-level numeric input as second argument.\n") ;
  }
  if (!mxIsInt32(prhs[1])) {
    mexErrMsgTxt("ml_haralicktexture requires as second input Gray Level vector of type int32 .\n") ;
  }
  if (!mxIsNumeric(prhs[2])) {
    mexErrMsgTxt("ml_haralicktexture requires a single numeric input as a third argument.\n") ;
  }


  NDims = mxGetNumberOfDimensions(prhs[0]);
  dims = mxGetDimensions(prhs[0]);
  if (NDims != 2) {
    mexErrMsgTxt("ml_3Dtexture requires a 2D matrix as the first input.\n");
    if (dims[0]!=dims[1]){
    	
    }    
  }
  coomat = (float*)mxGetData(prhs[0]) ;
  grayLevels = (int*)mxGetData(prhs[1]) ;
  nbGrayLevels = (int*)mxGetData(prhs[2]) ;
  /* mexPrintf("coomat=%f, %f, %f, %f, %f, %f, %f, %f, %f \n",*coomat,*(coomat+1),*(coomat+2),*(coomat+3),*(coomat+4),*(coomat+5),*(coomat+6),*(coomat+7),*(coomat+8));
     mexPrintf("grayLevels=%i, %i, %i \n",*grayLevels,*(grayLevels+1),*(grayLevels+2));
     mexPrintf("nbGrayLevels = %i\n",*nbGrayLevels); 
  */   

  features_used = mxCalloc(1, sizeof(TEXTURE_FEATURE_MAP)) ;
  if(!features_used) 
    mexErrMsgTxt("ml_haralicktexture: error allocating features_used.") ;

  features_used->ASM 	 	 = 1 ;
  features_used->contrast 	 = 1 ;
  features_used->correlation = 1 ;
  features_used->variance 	 = 1 ;
  features_used->IDM 		 = 1 ;
  features_used->sum_avg 	 = 1 ;
  features_used->sum_var	 = 1 ;
  features_used->sum_entropy = 1 ;
  features_used->entropy 	 = 1 ;
  features_used->diff_var 	 = 1 ;
  features_used->diff_entropy= 1 ;
  features_used->meas_corr1  = 1 ;
  features_used->meas_corr2  = 1 ;
  features_used->max_corr_coef = 1 ;
  
  features = mxCalloc(1, sizeof(TEXTURE));
  if (!features) mexErrMsgTxt("error allocating features.\n");
  
  feature_names = mxCalloc(14, sizeof(char*));
  for (i=0;i<15;i++)
  	feature_names[i] = mxCalloc(23, sizeof(char));

  Extract_Texture_Features(coomat,nbGrayLevels,grayLevels,features_used,features, feature_names); 
 
  outputsize[row] = 14 ;
  outputsize[col] = 1; 

  /* return the vector of texture features */
  plhs[0] = mxCreateNumericArray(2, outputsize, mxSINGLE_CLASS, mxREAL) ;
  if (!plhs[0]) mexErrMsgTxt("mb_texture: error allocating return variable.") ;
  
  output = (float*)mxGetData(plhs[0]) ;
  
  outputsize[row] = 14 ;
  outputsize[col] = 23; 
  
  plhs[1]= mxCreateCharMatrixFromStrings(outputsize[row], (const char **)feature_names); 


  /* Copy the features into the return variable */
	outputindex[col]=0 ;
    outputindex[row]=0 ;
    offset =  mxCalcSingleSubscript(plhs[0], 2, outputindex) ;
    output[offset] = features->ASM;

    outputindex[row]++ ;
    offset =  mxCalcSingleSubscript(plhs[0], 2, outputindex) ;
    output[offset] = features->contrast;

    outputindex[row]++ ;
    offset =  mxCalcSingleSubscript(plhs[0], 2, outputindex) ;
    output[offset] = features->correlation;

    outputindex[row]++ ;
    offset =  mxCalcSingleSubscript(plhs[0], 2, outputindex) ;
    output[offset] = features->variance;

    outputindex[row]++ ;
    offset =  mxCalcSingleSubscript(plhs[0], 2, outputindex) ;
    output[offset] = features->IDM;

    outputindex[row]++ ;
    offset =  mxCalcSingleSubscript(plhs[0], 2, outputindex) ;
    output[offset] = features->sum_avg;

    outputindex[row]++ ;
    offset =  mxCalcSingleSubscript(plhs[0], 2, outputindex) ;
    output[offset] = features->sum_var;

    outputindex[row]++ ;
    offset =  mxCalcSingleSubscript(plhs[0], 2, outputindex) ;
    output[offset] = features->sum_entropy;

    outputindex[row]++ ;
    offset =  mxCalcSingleSubscript(plhs[0], 2, outputindex) ;
    output[offset] = features->entropy;

    outputindex[row]++ ;
    offset =  mxCalcSingleSubscript(plhs[0], 2, outputindex) ;
    output[offset] = features->diff_var;

    outputindex[row]++ ;
    offset =  mxCalcSingleSubscript(plhs[0], 2, outputindex) ;
    output[offset] = features->diff_entropy;

    outputindex[row]++ ;
    offset =  mxCalcSingleSubscript(plhs[0], 2, outputindex) ;
    output[offset] = features->meas_corr1;

    outputindex[row]++ ;
    offset =  mxCalcSingleSubscript(plhs[0], 2, outputindex) ;
    output[offset] = features->meas_corr2;

    outputindex[row]++ ;
    offset =  mxCalcSingleSubscript(plhs[0], 2, outputindex) ;
    output[offset] = features->max_corr_coef;
    
  /*
    Memory clean-up.
  */
  mxFree(features_used) ;
  mxFree(features);
  for (i=0; i<14 ; i++)
    mxFree(feature_names[i]);
  mxFree(feature_names);
}
