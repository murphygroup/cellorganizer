/* ml_3Dpgmcoocmat.c - calculate the cooccurrence matrix
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
**  Nov 4, 2008 -  Estelle Glory : now this code calculate only the cooccurrence matrix and return it   
*/




int* Calculate_Cooccurrence_Matrix( register gray *grays, int nx, int ny, int nz, int** directions, int ndir, float** P_matrix, int *nbGraylevels,float* output, int *ab)  
{

  int tonec[PGM_MAXMAXVAL+1],R[13], x, y;
  int directionToCount[13];
  int argn, bps, padright, row, col, i, j, k, m, d;
  int itone, jtone, g_val, offset,tones;
  int *grayLevels;
  
  /*float **P_matrix; */
  gray nmaxval, maxval;

  d = 1; 





 

}

