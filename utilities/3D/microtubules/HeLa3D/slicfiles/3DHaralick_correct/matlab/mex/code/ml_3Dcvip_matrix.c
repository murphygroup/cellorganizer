/* ml_3Dpgmtxtur.c - calculate textural features on a portable graymap
**
** Author: James Darrell McCauley
**         Texas Agricultural Experiment Station
**         Department of Agricultural Engineering
**         Texas A&M University
**         College Station, Texas 77843-2117 USA
**
** Code written partially taken from pgmtofs.c in the PBMPLUS package
** by Jef Poskanzer.
**
** Algorithms for calculating features (and some explanatory comments) are
** taken from:
**
**   Haralick, R.M., K. Shanmugam, and I. Dinstein. 1973. Textural features
**   for image classification.  IEEE Transactions on Systems, Man, and
**   Cybertinetics, SMC-3(6):610-621.
**
** Copyright (C) 1991 Texas Agricultural Experiment Station, employer for
** hire of James Darrell McCauley
**
** Permission to use, copy, modify, and distribute this software and its
** documentation for any purpose and without fee is hereby granted, provided
** that the above copyright notice appear in all copies and that both that
** copyright notice and this permission notice appear in supporting
** documentation.  This software is provided "as is" without express or
** implied warranty.
**
** THE TEXAS AGRICULTURAL EXPERIMENT STATION (TAES) AND THE TEXAS A&M
** UNIVERSITY SYSTEM (TAMUS) MAKE NO EXPRESS OR IMPLIED WARRANTIES
** (INCLUDING BY WAY OF EXAMPLE, MERCHANTABILITY) WITH RESPECT TO ANY
** ITEM, AND SHALL NOT BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL
** OR CONSEQUENTAL DAMAGES ARISING OUT OF THE POSESSION OR USE OF
** ANY SUCH ITEM. LICENSEE AND/OR USER AGREES TO INDEMNIFY AND HOLD
** TAES AND TAMUS HARMLESS FROM ANY CLAIMS ARISING OUT OF THE USE OR
** POSSESSION OF SUCH ITEMS.
** 
** Modification History:
** 24 Jun 91 - J. Michael Carstensen <jmc@imsor.dth.dk> supplied fix for 
**             correlation function.
**
** Aug. 7 96 - Wenxing Li: huge memory leaks are fixed.
**
** 23 Nov 98 - M. Boland : Compile with the following for use with Matlab
**               under Red Hat Linux 5.1 (i.e. use libc5 instead of glibc)
**
**gcc -c -B /usr/libc5/usr/lib/gcc-lib/ -nostdinc -nostdinc++ -I/usr/libc5/usr/include -I/usr/libc5/usr/lib/gcc-lib/i386-linux/2.7.2.1/include -I/home/boland/**Matlab/Mex/Include -ansi cvip_pgmtexture.c
**
**ar r libmb_cvip.a cvip_pgmtexture.o
**
**>> mex -f gccopts.sh -lmb_cvip -L/home/boland/Matlab/Mex  mb_texture.c
**
**
**   29 Nov 98 - M. Boland : Modified calculations to produce the same values
**                 as the kharalick routine in Khoros.  Some feature 
**                 calculations below were wrong, others made different 
**                 assumptions (the Haralick paper is not always explicit).


    Aug 19, 2002 - Xiang Chen : Modified program to work on 3D images
    
**  November 1, 2008 - Estelle just calculate the cooccurrence matrix    
*/

#include <math.h>
#include "Include/ppgm.h"
#include "Include/3DCVIPtexture.h"
 
#define RADIX 2.0
#define EPSILON 0.000000001
#define BL  "Angle                 "
#define F1  "Angular Second Moment \0"
#define F2  "Contrast              \0"
#define F3  "Correlation           \0"
#define F4  "Variance              \0"
#define F5  "Inverse Diff Moment   \0"
#define F6  "Sum Average           \0"
#define F7  "Sum Variance          \0"
#define F8  "Sum Entropy           \0"
#define F9  "Entropy               \0"
#define F10 "Difference Variance   \0"
#define F11 "Difference Entropy    \0"
#define F12 "Meas of Correlation-1 \0"
#define F13 "Meas of Correlation-2 \0"
#define F14 "Max Correlation Coeff \0"
#define SIGN(x,y) ((y)<0 ? -fabs(x) : fabs(x))
#define DOT fprintf(stderr,".")
#define SWAP(a,b) {y=(a);(a)=(b);(b)=y;}
#define idx(x, y, z) (y) + (x) * ny + (z) * ny * nx

void results (),results2(),  mkbalanced (), reduction (), simplesrt ();
int hessenberg ();
float f1_asm (), f2_contrast (), f3_corr (), f4_var (), f5_idm (),
 f6_savg (), f7_svar (), f8_sentropy (), f9_entropy (), f10_dvar (),
 f11_dentropy (), f12_icorr (), f13_icorr (), f14_maxcorr (), *pgm_vector (),
 **pgm_matrix ();


void Extract_Texture_Features(float* coomat, int* nbGrayLevels, int* grayLevels, TEXTURE_FEATURE_MAP *feature_usage, TEXTURE *Texture,char **feat_name)
{
  FILE *ifp;
  register gray  *gP;
  int tonec[PGM_MAXMAXVAL+1], tone[PGM_MAXMAXVAL+1], angle, d = 1, x, y;
  int argn, bps, padright, row, col, i, j, k;
  float ASMTmp, contrastTmp, corrTmp, varTmp, idmTmp, savgTmp;
  float sentropyTmp, svarTmp, entropyTmp, dvarTmp, dentropyTmp;
  float icorrTmp, maxcorrTmp;
  int itone, jtone, g_val;
  float **P_matrix;
  float *Tp;
  gray nmaxval, maxval;
  char *usage = "[-d <d>] [pgmfile]";
 /* char **feat_name;*/

  strcpy(feat_name[0],F1);
  strcpy(feat_name[1],F2);
  strcpy(feat_name[2],F3);
  strcpy(feat_name[3],F4);
  strcpy(feat_name[4],F5);
  strcpy(feat_name[5],F6);
  strcpy(feat_name[6],F7);
  strcpy(feat_name[7],F8);
  strcpy(feat_name[8],F9);
  strcpy(feat_name[9],F10);
  strcpy(feat_name[10],F11);
  strcpy(feat_name[11],F12);
  strcpy(feat_name[12],F13);
  strcpy(feat_name[13],F14);
  
 
  /*tonec= grayLevels
    fprintf (stderr, " in 1 -coomat=%f, %f, %f, %f, %f, %f, %f, %f, %f \n",*coomat,*(coomat+1),*(coomat+2),*(coomat+3),*(coomat+4),*(coomat+5),*(coomat+6),*(coomat+7),*(coomat+8));
    fprintf (stderr, " in 1 -grayLevels=%i, %i, %i \n",*grayLevels,*(grayLevels+1),*(grayLevels+2));
    fprintf (stderr, " in 1 -nbGrayLevels = %i\n",*nbGrayLevels);
  */
  
  P_matrix= pgm_matrix (0, (*nbGrayLevels)-1, 0, (*nbGrayLevels)-1);
  for (j=0;j<*nbGrayLevels;j++)
    for (i=0;i<*nbGrayLevels;i++)
      P_matrix[j][i]=coomat[j*(*nbGrayLevels)+i];

  if (feature_usage->ASM)
  	ASMTmp = f1_asm (P_matrix, *nbGrayLevels);
  Tp = &Texture->ASM;
  results2 (Tp, F1, ASMTmp);

  if (feature_usage->contrast)
  	contrastTmp = f2_contrast (P_matrix, *nbGrayLevels);
  Tp = &Texture->contrast;
  results2 (Tp, F2, contrastTmp);

  if (feature_usage->correlation)
  	corrTmp = f3_corr (P_matrix, *nbGrayLevels);
  Tp = &Texture->correlation;
  results2 (Tp, F3, corrTmp);

  if (feature_usage->variance)
  	varTmp = f4_var (P_matrix, *nbGrayLevels);
  Tp = &Texture->variance;
  results2 (Tp, F4, varTmp); 

  if (feature_usage->IDM)
  	idmTmp = f5_idm (P_matrix, *nbGrayLevels);
  Tp = &Texture->IDM;
  results2 (Tp, F5, idmTmp); 

  if (feature_usage->sum_avg)
  	savgTmp = f6_savg (P_matrix, *nbGrayLevels);
  Tp = &Texture->sum_avg;
  results2 (Tp, F6, savgTmp); 

  if (feature_usage->sum_var)
  	sentropyTmp = f8_sentropy (P_matrix, *nbGrayLevels);
  if (feature_usage->sum_entropy)  
  	svarTmp = f7_svar (P_matrix, *nbGrayLevels, sentropyTmp);
  Tp = &Texture->sum_var;
  results2 (Tp, F7, svarTmp); 
  Tp = &Texture->sum_entropy;
  results2 (Tp, F8, sentropyTmp); 

  if (feature_usage->entropy)
  	entropyTmp = f9_entropy (P_matrix, *nbGrayLevels);
  Tp = &Texture->entropy;
  results2 (Tp, F9, entropyTmp); 

  if (feature_usage->diff_var)
	dvarTmp = f10_dvar (P_matrix, *nbGrayLevels);
  Tp = &Texture->diff_var;
  results2 (Tp, F10, dvarTmp);

  if (feature_usage->diff_entropy)
	dentropyTmp = f11_dentropy (P_matrix, *nbGrayLevels);
  Tp = &Texture->diff_entropy;
  results2 (Tp, F11, dentropyTmp);

   if (feature_usage->meas_corr1)
	icorrTmp = f12_icorr (P_matrix, *nbGrayLevels);
  Tp = &Texture->meas_corr1;
  results2 (Tp, F12, icorrTmp);

  if (feature_usage->meas_corr2)
	icorrTmp = f13_icorr (P_matrix, *nbGrayLevels);
  Tp = &Texture->meas_corr2;
  results2 (Tp, F13, icorrTmp);

  if (feature_usage->max_corr_coef) 
	maxcorrTmp = f14_maxcorr (P_matrix, *nbGrayLevels);
  /* M. Boland - 24 Nov 98 */
  else 
    maxcorrTmp = 0 ;

  Tp = &Texture->max_corr_coef;
  results2 (Tp, F14, maxcorrTmp);

  for (j=0; j<(*nbGrayLevels); j++)
    free(P_matrix[j]);

  return;
}

float f1_asm (P, Ng)
  float **P;
  int Ng;

/* Angular Second Moment */
{
  int i, j;
  float sum = 0;

  for (i = 0; i < Ng; ++i)
    for (j = 0; j < Ng; ++j)
      sum += P[i][j] * P[i][j];

  return sum;

  /*
   * The angular second-moment feature (ASM) f1 is a measure of homogeneity
   * of the image. In a homogeneous image, there are very few dominant
   * gray-tone transitions. Hence the P matrix for such an image will have
   * fewer entries of large magnitude.
   */
}


float f2_contrast (P, Ng)
  float **P;
  int Ng;

/* Contrast */
{
  int i, j, n;
  float sum = 0, bigsum = 0;

  for (n = 0; n < Ng; ++n)
  {
    for (i = 0; i < Ng; ++i)
      for (j = 0; j < Ng; ++j)
	if ((i - j) == n || (j - i) == n)
	  sum += P[i][j];
    bigsum += n * n * sum;

    sum = 0;
  }
  return bigsum;

  /*
   * The contrast feature is a difference moment of the P matrix and is a
   * measure of the contrast or the amount of local variations present in an
   * image.
   */
}

float f3_corr (P, Ng)
  float **P;
  int Ng;

/* Correlation */
{
  int i, j;
  float sum_sqrx = 0, sum_sqry = 0, tmp, *px;
  float meanx =0 , meany = 0 , stddevx, stddevy;

  px = pgm_vector (0, Ng);
  for (i = 0; i < Ng; ++i)
    px[i] = 0;

  /*
   * px[i] is the (i-1)th entry in the marginal probability matrix obtained
   * by summing the rows of p[i][j]
   */
  for (i = 0; i < Ng; ++i)
    for (j = 0; j < Ng; ++j)
      px[i] += P[i][j];


  /* Now calculate the means and standard deviations of px and py */
  /*- fix supplied by J. Michael Christensen, 21 Jun 1991 */
  /*- further modified by James Darrell McCauley, 16 Aug 1991 
   *     after realizing that meanx=meany and stddevx=stddevy
   */
  for (i = 0; i < Ng; ++i)
  {
    meanx += px[i]*i;
    sum_sqrx += px[i]*i*i;
  }
  /* M. Boland meanx = meanx/(sqrt(Ng)); */
  meany = meanx;
  sum_sqry = sum_sqrx;
  stddevx = sqrt (sum_sqrx - (meanx * meanx));
  stddevy = stddevx;

  /* Finally, the correlation ... */
  for (tmp = 0, i = 0; i < Ng; ++i)
    for (j = 0; j < Ng; ++j)
      tmp += i*j*P[i][j];

 free(px); 
 return (tmp - meanx * meany) / (stddevx * stddevy);
  /*
   * This correlation feature is a measure of gray-tone linear-dependencies
   * in the image.
   */
}


float f4_var (P, Ng)
  float **P;
  int Ng;

/* Sum of Squares: Variance */
{
  int i, j;
  float mean = 0, var = 0;

  /*- Corrected by James Darrell McCauley, 16 Aug 1991
   *  calculates the mean intensity level instead of the mean of
   *  cooccurrence matrix elements 
   */
  for (i = 0; i < Ng; ++i)
    for (j = 0; j < Ng; ++j)
      mean += i * P[i][j];

  for (i = 0; i < Ng; ++i)
    for (j = 0; j < Ng; ++j)
      /*  M. Boland - var += (i + 1 - mean) * (i + 1 - mean) * P[i][j]; */
      var += (i - mean) * (i - mean) * P[i][j];

  return var;
}

float f5_idm (P, Ng)
  float **P;
  int Ng;

/* Inverse Difference Moment */
{
  int i, j;
  float idm = 0;

  for (i = 0; i < Ng; ++i)
    for (j = 0; j < Ng; ++j)
      idm += P[i][j] / (1 + (i - j) * (i - j));

  return idm;
}

float Pxpy[2 * PGM_MAXMAXVAL];

float f6_savg (P, Ng)
  float **P;
  int Ng;

/* Sum Average */
{
  int i, j;
  extern float Pxpy[2 * PGM_MAXMAXVAL];
  float savg = 0;

  for (i = 0; i <= 2 * Ng; ++i)
    Pxpy[i] = 0;

  for (i = 0; i < Ng; ++i)
    for (j = 0; j < Ng; ++j)
      /* M. Boland Pxpy[i + j + 2] += P[i][j]; */
      /* Indexing from 2 instead of 0 is inconsistent with rest of code*/
      Pxpy[i + j] += P[i][j];
  /* M. Boland for (i = 2; i <= 2 * Ng; ++i) */
  /* Indexing from 2 instead of 0 is inconsistent with rest of code*/
  for (i = 0; i <= (2 * Ng - 2); ++i)
    savg += i * Pxpy[i];

  return savg;
}


float f7_svar (P, Ng, S)
  float **P, S;
  int Ng;

/* Sum Variance */
{
  int i, j;
  extern float Pxpy[2 * PGM_MAXMAXVAL];
  float var = 0;

  for (i = 0; i <= 2 * Ng; ++i)
    Pxpy[i] = 0;

  for (i = 0; i < Ng; ++i)
    for (j = 0; j < Ng; ++j)
      /* M. Boland Pxpy[i + j + 2] += P[i][j]; */
      /* Indexing from 2 instead of 0 is inconsistent with rest of code*/
      Pxpy[i + j] += P[i][j];

  /*  M. Boland for (i = 2; i <= 2 * Ng; ++i) */
  /* Indexing from 2 instead of 0 is inconsistent with rest of code*/
  for (i = 0; i <= (2 * Ng - 2); ++i)
    var += (i - S) * (i - S) * Pxpy[i];

  return var;
}

float f8_sentropy (P, Ng)
  float **P;
  int Ng;

/* Sum Entropy */
{
  int i, j;
  extern float Pxpy[2 * PGM_MAXMAXVAL];
  float sentropy = 0;

  for (i = 0; i <= 2 * Ng; ++i)
    Pxpy[i] = 0;

  for (i = 0; i < Ng; ++i)
    for (j = 0; j < Ng; ++j)
      Pxpy[i + j + 2] += P[i][j];

  for (i = 2; i <= 2 * Ng; ++i)
    /*  M. Boland  sentropy -= Pxpy[i] * log10 (Pxpy[i] + EPSILON); */
    sentropy -= Pxpy[i] * log10 (Pxpy[i] + EPSILON)/log10(2.0) ;

  return sentropy;
}


float f9_entropy (P, Ng)
  float **P;
  int Ng;

/* Entropy */
{
  int i, j;
  float entropy = 0;

  for (i = 0; i < Ng; ++i)
    for (j = 0; j < Ng; ++j)
      /*      entropy += P[i][j] * log10 (P[i][j] + EPSILON); */
      entropy += P[i][j] * log10 (P[i][j] + EPSILON)/log10(2.0) ;

  return -entropy; 
}


float f10_dvar (P, Ng)
  float **P;
  int Ng;

/* Difference Variance */
{
  int i, j, tmp;
  extern float Pxpy[2 * PGM_MAXMAXVAL];
  float sum = 0, sum_sqr = 0, var = 0;

  for (i = 0; i <= 2 * Ng; ++i)
    Pxpy[i] = 0;

  for (i = 0; i < Ng; ++i)
    for (j = 0; j < Ng; ++j)
      Pxpy[abs (i - j)] += P[i][j];

  /* Now calculate the variance of Pxpy (Px-y) */
  for (i = 0; i < Ng; ++i)
  {
    sum += i * Pxpy[i] ;
    sum_sqr += i * i * Pxpy[i] ;
    /* M. Boland sum += Pxpy[i];
    sum_sqr += Pxpy[i] * Pxpy[i];*/
  }
  /*tmp = Ng * Ng ;  M. Boland - wrong anyway, should be Ng */
  /*var = ((tmp * sum_sqr) - (sum * sum)) / (tmp * tmp); */
  
  var = sum_sqr - sum*sum ;

  return var;
}

float f11_dentropy (P, Ng)
  float **P;
  int Ng;

/* Difference Entropy */
{
  int i, j, tmp;
  extern float Pxpy[2 * PGM_MAXMAXVAL];
  float sum = 0, sum_sqr = 0, var = 0;

  for (i = 0; i <= 2 * Ng; ++i)
    Pxpy[i] = 0;

  for (i = 0; i < Ng; ++i)
    for (j = 0; j < Ng; ++j)
      Pxpy[abs (i - j)] += P[i][j];

  for (i = 0; i < Ng; ++i)
    /*    sum += Pxpy[i] * log10 (Pxpy[i] + EPSILON); */
    sum += Pxpy[i] * log10 (Pxpy[i] + EPSILON)/log10(2.0) ;

  return -sum;
}

float f12_icorr (P, Ng)
  float **P;
  int Ng;

/* Information Measures of Correlation */
/* All /log10(2.0) added by M. Boland */
{
  int i, j, tmp;
  float *px, *py;
  float hx = 0, hy = 0, hxy = 0, hxy1 = 0, hxy2 = 0;

  px = pgm_vector (0, Ng);
  py = pgm_vector (0, Ng);

  /*
   * px[i] is the (i-1)th entry in the marginal probability matrix obtained
   * by summing the rows of p[i][j]
   */
  for (i = 0; i < Ng; ++i)
  {
    for (j = 0; j < Ng; ++j)
    {
      px[i] += P[i][j];
      py[j] += P[i][j];
    }
  }

  for (i = 0; i < Ng; ++i)
    for (j = 0; j < Ng; ++j)
    {
      hxy1 -= P[i][j] * log10 (px[i] * py[j] + EPSILON)/log10(2.0);
      hxy2 -= px[i] * py[j] * log10 (px[i] * py[j] + EPSILON)/log10(2.0);
      hxy -= P[i][j] * log10 (P[i][j] + EPSILON)/log10(2.0);
    }

  /* Calculate entropies of px and py - is this right? */
  for (i = 0; i < Ng; ++i)
  {
    hx -= px[i] * log10 (px[i] + EPSILON)/log10(2.0);
    hy -= py[i] * log10 (py[i] + EPSILON)/log10(2.0);
  }
/*  fprintf(stderr,"hxy1=%f\thxy=%f\thx=%f\thy=%f\n",hxy1,hxy,hx,hy); */
  free(px);
  free(py);
  return ((hxy - hxy1) / (hx > hy ? hx : hy));
}

float f13_icorr (P, Ng)
  float **P;
  int Ng;

/* Information Measures of Correlation */
/* All /log10(2.0) added by M. Boland */

{
  int i, j;
  float *px, *py;
  float hx = 0, hy = 0, hxy = 0, hxy1 = 0, hxy2 = 0;

  px = pgm_vector (0, Ng);
  py = pgm_vector (0, Ng);

  /*
   * px[i] is the (i-1)th entry in the marginal probability matrix obtained
   * by summing the rows of p[i][j]
   */
  for (i = 0; i < Ng; ++i)
  {
    for (j = 0; j < Ng; ++j)
    {
      px[i] += P[i][j];
      py[j] += P[i][j];
    }
  }

  for (i = 0; i < Ng; ++i)
    for (j = 0; j < Ng; ++j)
    {
      hxy1 -= P[i][j] * log10 (px[i] * py[j] + EPSILON)/log10(2.0);
      hxy2 -= px[i] * py[j] * log10 (px[i] * py[j] + EPSILON)/log10(2.0);
      hxy -= P[i][j] * log10 (P[i][j] + EPSILON)/log10(2.0);
    }

  /* Calculate entropies of px and py */
  for (i = 0; i < Ng; ++i)
  {
    hx -= px[i] * log10 (px[i] + EPSILON)/log10(2.0);
    hy -= py[i] * log10 (py[i] + EPSILON)/log10(2.0);
  }
/*  fprintf(stderr,"hx=%f\thxy2=%f\n",hx,hxy2); */
  free(px);
  free(py);
  return (sqrt (abs (1 - exp (-2.0 * (hxy2 - hxy)))));
}

float f14_maxcorr (P, Ng)
  float **P;
  int Ng;

/* Returns the Maximal Correlation Coefficient */
{
  int i, j, k;
  float *px, *py, **Q;
  float *x, *iy, tmp;
  float f;

  px = pgm_vector (0, Ng);
  py = pgm_vector (0, Ng);
  Q = pgm_matrix (1, Ng + 1, 1, Ng + 1);
  x = pgm_vector (1, Ng);
  iy = pgm_vector (1, Ng);

  /*
   * px[i] is the (i-1)th entry in the marginal probability matrix obtained
   * by summing the rows of p[i][j]
   */
  for (i = 0; i < Ng; ++i)
  {
    for (j = 0; j < Ng; ++j)
    {
      px[i] += P[i][j];
      py[j] += P[i][j];
    }
  }

  /* Find the Q matrix */
  for (i = 0; i < Ng; ++i)
  {
    for (j = 0; j < Ng; ++j)
    {
      Q[i + 1][j + 1] = 0;
      for (k = 0; k < Ng; ++k)
	Q[i + 1][j + 1] += P[i][k] * P[j][k] / px[i] / py[k];
    }
  }

  /* Balance the matrix */
  mkbalanced (Q, Ng);
  /* Reduction to Hessenberg Form */
  reduction (Q, Ng);
  /* Finding eigenvalue for nonsymetric matrix using QR algorithm */
  if (!hessenberg (Q, Ng, x, iy))
	{ for (i=1; i<=Ng+1; i++) free(Q[i]+1);
	  free(Q+1);
	  free((char *)px);
	  free((char *)py);
	  free((x+1));
	  free((iy+1));
	  return 0.0;
	  /* fixed for Linux porting,
	   * I don't know what should be returned
	   */
	}
  /* simplesrt(Ng,x); */
  /* Returns the sqrt of the second largest eigenvalue of Q */
  for (i = 2, tmp = x[1]; i <= Ng; ++i)
    tmp = (tmp > x[i]) ? tmp : x[i];

  f = sqrt(x[Ng - 1]);

 for (i=1; i<=Ng+1; i++) free(Q[i]+1);
 free(Q+1);
 free((char *)px); 
 free((char *)py); 
 free((x+1)); 
 free((iy+1)); 

 return f;
}

float *pgm_vector (nl, nh)
  int nl, nh;
{
  float *v;
  int    i;

  v = (float *) malloc ((unsigned) (nh - nl + 1) * sizeof (float));
  if (!v)
    fprintf (stderr, "memory allocation failure (pgm_vector) "), exit (1);

  for (i=0; i<=(nh-nl); i++) v[i]=0;
  return v - nl;
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

void results (Tp, c, a)
  float *Tp;
  char *c;
  float *a;
{
  int i;
  float max, min, sum;
  max = a[0];
  min = a[0];
  sum = 0;
/*  DOT;
  fprintf (stdout, "%s", c);
*/  for (i = 0; i < 13; ++i, *Tp++)
    {	
    if (a[i] <= min)
	min = a[i];
    if (a[i] > max)
	max = a[i];
  /*  fprintf (stdout, "% 1.3e ", a[i]); */
    *Tp = a[i];
    sum += a[i];
    }	
/*  fprintf (stdout, "% 1.3e  % 1.3e\n", (a[0] + a[1] + a[2] + a[3]) / 4,max-min); */
  *Tp = sum / 13;
  *Tp++;
  *Tp = max - min;
 
  	
}
void results2 (Tp, c, a)
  float *Tp;
  char *c;
  float a;
{
  int i;
 
/*  DOT;
  fprintf (stdout, "%s", c);
*/  
  /*  fprintf (stdout, "% 1.3e ", a[i]); */
   *Tp = a;
/*  fprintf (stdout, "% 1.3e  % 1.3e\n", (a[0] + a[1] + a[2] + a[3]) / 4,max-min); */
 	
}

void simplesrt (n, arr)
  int n;
  float arr[];
{
  int i, j;
  float a;

  for (j = 2; j <= n; j++)
  {
    a = arr[j];
    i = j - 1;
    while (i > 0 && arr[i] > a)
    {
      arr[i + 1] = arr[i];
      i--;
    }
    arr[i + 1] = a;
  }
}

void mkbalanced (a, n)
  float **a;
  int n;
{
  int last, j, i;
  float s, r, g, f, c, sqrdx;

  sqrdx = RADIX * RADIX;
  last = 0;
  while (last == 0)
  {
    last = 1;
    for (i = 1; i <= n; i++)
    {
      r = c = 0.0;
      for (j = 1; j <= n; j++)
	if (j != i)
	{
	  c += fabs (a[j][i]);
	  r += fabs (a[i][j]);
	}
      if (c && r)
      {
	g = r / RADIX;
	f = 1.0;
	s = c + r;
	while (c < g)
	{
	  f *= RADIX;
	  c *= sqrdx;
	}
	g = r * RADIX;
	while (c > g)
	{
	  f /= RADIX;
	  c /= sqrdx;
	}
	if ((c + r) / f < 0.95 * s)
	{
	  last = 0;
	  g = 1.0 / f;
	  for (j = 1; j <= n; j++)
	    a[i][j] *= g;
	  for (j = 1; j <= n; j++)
	    a[j][i] *= f;
	}
      }
    }
  }
}


void reduction (a, n)
  float **a;
  int n;
{
  int m, j, i;
  float y, x;

  for (m = 2; m < n; m++)
  {
    x = 0.0;
    i = m;
    for (j = m; j <= n; j++)
    {
      if (fabs (a[j][m - 1]) > fabs (x))
      {
	x = a[j][m - 1];
	i = j;
      }
    }
    if (i != m)
    {
      for (j = m - 1; j <= n; j++)
	SWAP (a[i][j], a[m][j])  
	for (j = 1; j <= n; j++)
	  SWAP (a[j][i], a[j][m]) 
	  a[j][i] = a[j][i];
    }
    if (x)
    {
      for (i = m + 1; i <= n; i++)
      {
	if (y = a[i][m - 1])
	{
	  y /= x;
	  a[i][m - 1] = y;
	  for (j = m; j <= n; j++)
	    a[i][j] -= y * a[m][j];
	  for (j = 1; j <= n; j++)
	    a[j][m] += y * a[j][i];
	}
      }
    }
  }
}

int hessenberg (a, n, wr, wi)
  float **a, wr[], wi[];
  int n;

{
  int nn, m, l, k, j, its, i, mmin;
  float z, y, x, w, v, u, t, s, r, q, p, anorm;

  anorm = fabs (a[1][1]);
  for (i = 2; i <= n; i++)
    for (j = (i - 1); j <= n; j++)
      anorm += fabs (a[i][j]);
  nn = n;
  t = 0.0;
  while (nn >= 1)
  {
    its = 0;
    do
    {
      for (l = nn; l >= 2; l--)
      {
	s = fabs (a[l - 1][l - 1]) + fabs (a[l][l]);
	if (s == 0.0)
	  s = anorm;
	if ((float) (fabs (a[l][l - 1]) + s) == s)
	  break;
      }
      x = a[nn][nn];
      if (l == nn)
      {
	wr[nn] = x + t;
	wi[nn--] = 0.0;
      }
      else
      {
	y = a[nn - 1][nn - 1];
	w = a[nn][nn - 1] * a[nn - 1][nn];
	if (l == (nn - 1))
	{
	  p = 0.5 * (y - x);
	  q = p * p + w;
	  z = sqrt (fabs (q));
	  x += t;
	  if (q >= 0.0)
	  {
	    z = p + SIGN (z, p); 
	    wr[nn - 1] = wr[nn] = x + z;
	    if (z)
	      wr[nn] = x - w / z;
	    wi[nn - 1] = wi[nn] = 0.0;
	  }
	  else
	  {
	    wr[nn - 1] = wr[nn] = x + p;
	    wi[nn - 1] = -(wi[nn] = z);
	  }
	  nn -= 2;
	}
	else
	{
	  if (its == 30)
	    {
/*	    fprintf (stderr, 
                    "Too many iterations to required to find %s\ngiving up\n", 
                     F14);  */
	     return 0; /*exit (1);*/
	     }			
	  if (its == 10 || its == 20)
	  {
	    t += x;
	    for (i = 1; i <= nn; i++)
	      a[i][i] -= x;
	    s = fabs (a[nn][nn - 1]) + fabs (a[nn - 1][nn - 2]);
	    y = x = 0.75 * s;
	    w = -0.4375 * s * s;
	  }
	  ++its;
	  for (m = (nn - 2); m >= l; m--)
	  {
	    z = a[m][m];
	    r = x - z;
	    s = y - z;
	    p = (r * s - w) / a[m + 1][m] + a[m][m + 1];
	    q = a[m + 1][m + 1] - z - r - s;
	    r = a[m + 2][m + 1];
	    s = fabs (p) + fabs (q) + fabs (r);
	    p /= s;
	    q /= s;
	    r /= s;
	    if (m == l)
	      break;
	    u = fabs (a[m][m - 1]) * (fabs (q) + fabs (r));
	    v = fabs (p) * (fabs (a[m - 1][m - 1]) + fabs (z) + fabs (a[m + 1][m + 1]));
	    if ((float) (u + v) == v)
	      break;
	  }
	  for (i = m + 2; i <= nn; i++)
	  {
	    a[i][i - 2] = 0.0;
	    if (i != (m + 2))
	      a[i][i - 3] = 0.0;
	  }
	  for (k = m; k <= nn - 1; k++)
	  {
	    if (k != m)
	    {
	      p = a[k][k - 1];
	      q = a[k + 1][k - 1];
	      r = 0.0;
	      if (k != (nn - 1))
		r = a[k + 2][k - 1];
	      if (x = fabs (p) + fabs (q) + fabs (r))
	      {
		p /= x;
		q /= x;
		r /= x;
	      }
	    }
	    if (s = SIGN (sqrt (p * p + q * q + r * r), p)) 
	    {
	      if (k == m)
	      {
		if (l != m)
		  a[k][k - 1] = -a[k][k - 1];
	      }
	      else
		a[k][k - 1] = -s * x;
	      p += s;
	      x = p / s;
	      y = q / s;
	      z = r / s;
	      q /= p;
	      r /= p;
	      for (j = k; j <= nn; j++)
	      {
		p = a[k][j] + q * a[k + 1][j];
		if (k != (nn - 1))
		{
		  p += r * a[k + 2][j];
		  a[k + 2][j] -= p * z;
		}
		a[k + 1][j] -= p * y;
		a[k][j] -= p * x;
	      }
	      mmin = nn < k + 3 ? nn : k + 3;
	      for (i = l; i <= mmin; i++)
	      {
		p = x * a[i][k] + y * a[i][k + 1];
		if (k != (nn - 1))
		{
		  p += z * a[i][k + 2];
		  a[i][k + 2] -= p * r;
		}
		a[i][k + 1] -= p * q;
		a[i][k] -= p;
	      }
	    }
	  }
	}
      }
    } while (l < nn - 1);
  }
return 1;
}
