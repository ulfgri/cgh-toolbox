/*
 * A mex function for fast evaluation of polynomials
 * 
 * This software was developed at the National Institute of Standards and Technology
 * by employees of the Federal Government in the course of their official duties.
 * Pursuant to title 17 Section 105 of the United States Code this software is not
 * subject to copyright protection and is in the public domain. This software is an
 * experimental system. NIST assumes no responsibility whatsoever for its use by other
 * parties, and makes no guarantees, expressed or implied, about its quality, reliability,
 * or any other characteristic.
 *  
 * Short version: This software is in the PUBLIC DOMAIN.
 * Ulf Griesmann, June 2013, in better vectorizable form: July 2013
 */

#include <math.h>
#include <string.h>
#include "mex.h"

#define MAX(a,b) (((a)>(b))?(a):(b))

#ifdef __GNUC__
   #define RESTRICT __restrict
   #define INLINE __inline__
#else
   #define RESTRICT
   #define INLINE
#endif


/*-- local prototypes -----------------------------------------*/

static INLINE void
multiply_add(int N, double *a, double *b, double c);


/*-------------------------------------------------------------*/

void 
mexFunction(int nlhs, mxArray *plhs[], 
	    int nrhs, const mxArray *prhs[])
{
   int Ma, Na, Mx, Nx, La, Lx, k;
   double *A, *X, *Y;

   /* 
    * check number of arguments
    */
   if (nrhs != 2) {
      mexErrMsgTxt("fast_poly_eval :  exactly two arguments required.");
   }

   /*
    * check and get coefficient vector
    */
   Na = mxGetN(prhs[0]);
   Ma = mxGetM(prhs[0]);
   La = MAX(Na, Ma);  /* number of coefficients */
   A = mxGetData(prhs[0]);

   /*
    * check and get location vector
    */
   Nx = mxGetN(prhs[1]);
   Mx = mxGetM(prhs[1]);
   Lx = MAX(Nx, Mx);  /* number of evaluation points */
   X = mxGetData(prhs[1]);

   /* 
    * create output array
    */
   plhs[0] = mxCreateDoubleMatrix(Mx, Nx, mxREAL);
   Y = (double *)mxGetData(plhs[0]); /* Y is preset to 0 */

   /*
    * evaluate polynomial using Horner's algorithm
    */
   for (k=La-1; k>=0; k--)
      multiply_add(Lx, Y, X, A[k]);
}


/*
 * a fused multiply-add that the compiler can vectorize
 *
 *    a = a .* b + c
 *
 * a and b are vectors of length N, c is a scalar
 */
static INLINE void
multiply_add(int N, double * RESTRICT a, double * RESTRICT b, double c)
{
   int k;

   for (k=0; k<N; k++) {
      a[k] = a[k] * b[k] + c;
   }
}
