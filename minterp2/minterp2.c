/*
 * An Octave / Matlab mex function for very fast interpolation
 * of gridded data.
 *
 * For Octave, compile with:
 *    export CFLAGS="-O3 -march=native -fomit-frame-pointer" 
 *    mkoctfile --mex -s minterp2.c
 *
 * This software was developed at the National Institute of Standards and
 * Technology by employees of the United States Federal Government in the
 * course of their official duties. Pursuant to title 17 Section 105 of
 * the United States Code this software is not subject to copyright
 * protection and is in the public domain. This software is an
 * experimental system. NIST assumes no responsibility whatsoever for its
 * use by other parties, and makes no guarantees, expressed or implied,
 * about its quality, reliability, or any other characteristic.
 *
 * Ulf Griesmann, NIST, November 2011
 * add inline and restrict keywords to help compilers; U.G., July 2013
 */

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "mex.h"

#define STR_LEN    16
#define MAX(a, b) (((a)>(b)) ? (a) : (b))

#ifdef __GNUC__
   #define RESTRICT __restrict
   #define INLINE __inline__
#else
   #define RESTRICT
   #define INLINE
#endif

#define ERR_LEN 64


/*--- local function prototypes -------------------------------*/

static double interp_linear(double * RESTRICT pM, double x, double y);
static double interp_cubic(double * RESTRICT pM, double x, double y);
static double interp_pchip(double * RESTRICT pM, double x, double y);
static double interp_lagrange(double * RESTRICT pM, double x, double y);
static INLINE double dot_prod(int n, double * RESTRICT a, double * RESTRICT b);
static INLINE double dot_prodw(int n, double * RESTRICT a, 
			       double * RESTRICT b, double * RESTRICT w);

/* hash function */
#include "minterp2.h"


/*--- static variables ----------------------------------------*/

static double *map = NULL;  /* pointer to stored matrix (map) */
static double NaN;          /* system NaN representation */
static mwSize nr, nc, tnr;  /* rows, columns */
static mwSize mapsiz;       /* number of values in map */

static int int_init;
static mwIndex x0_last, y0_last;
static double ac[15];       /* polynomial coefficients */


/*-- Cleanup function ---------------------------------------------*/

void
cleanup_at_exit(void)
{
  if (map != NULL)
     mxFree(map);
}


/*-------------------------------------------------------------*/

void 
mexFunction(int nlhs, mxArray *plhs[], 
	    int nrhs, const mxArray *prhs[])
{
   mwSize np, npm, npn;
   int k;
   struct keyword *pk;  /* pointer to keyword structure */
   double *pM;          /* pointer to matrix (map) data */
   double *xi, *yi;     /* interpolation locations */   
   double *z;           /* vector with interpolation results */
   char imstr[STR_LEN]; /* interpolation method */
   char errmsg[ERR_LEN];
 

   /*************************************
    * store matrix if only one argument *
    *************************************/
   if (nrhs == 1) { 
      pM = mxGetData(prhs[0]); /* pointer to matrix data */
      nr = mxGetM(prhs[0]);    /* rows */
      tnr = 2*nr;
      nc = mxGetN(prhs[0]);    /* columns */

      /*
       * copy map to a persistent buffer
       */
      if (map == NULL) {          /* no buffer is allocated */
	 mapsiz = nr * nc;
	 map = mxCalloc(mapsiz, sizeof(double));
	 if (map == NULL)
	    mexErrMsgTxt("minterp2 :  failed to allocate memory for map.");
	 mexMakeMemoryPersistent(map);
	 mexAtExit(cleanup_at_exit);
      }
      else {                      /* a buffer already exists */
	 if (nr * nc != mapsiz) { /* but its size is wrong */
	    mapsiz = nr * nc;
	    map = mxRealloc(map, mapsiz*sizeof(double));
	    if (map == NULL)
	       mexErrMsgTxt("minterp2 :  failed to allocate memory for map.");
	    mexMakeMemoryPersistent(map);
	 }
      }

      /* 
       * Copy map data to the buffer.
       * The matrix is stored column-by-column (like Matlab)
       */
      memcpy(map, pM, mapsiz*sizeof(double));

      return;
   }

   /*
    * get system NaN
    */
   NaN = mxGetNaN();

   /*****************
    * interpolation *
    *****************/

   /* check if enough arguments were passed */
   if (nrhs < 3)
      mexErrMsgTxt("minterp2 :  expecting at least three input arguments.");

   /* check if first argument is empty */
   if ( mxIsEmpty(prhs[0]) ) { /* use the stored matrix */
      if (map == NULL)
	 mexErrMsgTxt("minterp2 :  no stored matrix.");
      pM = map;
   }
   else {                      /* use the argument matrix */
      if (map != NULL) {       /* free an existing stored map */
         mxFree(map);
	 map = NULL;
      }
      pM = mxGetData(prhs[0]);
      nr = mxGetM(prhs[0]);    /* rows */
      nc = mxGetN(prhs[0]);    /* columns */
      tnr = 2*nr;
   }
   
   xi = mxGetData(prhs[1]);    /* ptr to x values */
   yi = mxGetData(prhs[2]);    /* ptr to y values */
   npm = mxGetM(prhs[1]);      /* shape of xi */
   npn = mxGetN(prhs[1]);
   np = MAX(npm, npn);         /* number of points */

   /*
    * create output vector with the same shape as xi
    */
   plhs[0] = mxCreateDoubleMatrix(npm, npn, mxREAL);
   z = mxGetData(plhs[0]);

   /*
    * interpolate points
    */
   if (nrhs < 4) { /* interpolation set default method */
      strcpy(imstr, "pchip");
   }
   else {
      mxGetString(prhs[3], imstr, STR_LEN);
   }

   /*
    * look up the interpolation method
    */
   pk = (struct keyword *)in_word_set(imstr, strlen(imstr));
   if (pk == NULL) {
      sprintf(errmsg, "minterp2 :  unrecognized interpolation method -> %s", imstr);
      mexErrMsgTxt(errmsg);
   }

   int_init = 0;
   
   /*
    * dispatch selected interpolation function
    */
   for (k=0; k<np; k++)
      z[k] = (*pk->interp_func)(pM, xi[k], yi[k]);

   return;
}


/*-----------------------------------------------------------------*/

/*
 * bi-linear interpolation using the map values of four nodes surrounding 
 * the interpolation point
 */
static double 
interp_linear(double * RESTRICT pM, double x, double y) 
{
   register double z0, z1, z2, z3;
   register double X,Y;
   double xy[3];
   mwIndex x0, y0, idx00; 

   /* node (0,0) */
   x0 = floor(x);
   y0 = floor(y);
   
   /* return NaN if too close to the edge */
   if (x0 < 1 || x0 >= nc || y0 < 1 || y0 >= nr) {
      return NaN;
   }

   /* linear index of node (0,0) in map */
   idx00 = (x0-1)*nr + y0 -1;
   z0 = pM[idx00];

   /* evalute the polynomial at (x,y) */
   X = x - x0;
   Y = y - y0;
   if ( (X == 0.0) && (Y == 0.0) )
      return z0;

   /* re-calculate coefficients only if necessary */
   if (!int_init || x0 != x0_last || y0 != y0_last) {
   
      /* matrix values at nodes other than (0,0)*/
      z1 = pM[idx00 + nr];     /* node (1,0) */
      z2 = pM[idx00 + nr + 1]; /* node (1,1) */
      z3 = pM[idx00 + 1];      /* node (0,1) */

      /* solution for coefficients */
      ac[0] = z1 - z0;
      ac[1] = z3 - z0;
      ac[2] = z0 - z1 + z2 - z3;
      
      int_init = 1;
   }

   /* store node 0 location */
   x0_last = x0;
   y0_last = y0;

   xy[0] = X;
   xy[1] = Y;
   xy[2] = X*Y;

   return z0 + dot_prod(3,ac,xy);   
}


/*-----------------------------------------------------------------*/

/*
 * interpolate using the map values of 4 nodes surrounding
 * the interpolation point and the derivatives at the
 * nodes.
 */
static double 
interp_cubic(double * RESTRICT pM, double x, double y) 
{
   double z0, z1, z2, z3, z4, z5, z6, z7, z8, z9, z10, z11, z12, z13, z14, z15;
   double n4, n5, n6, n7, n8, n9, n10, n11, n12, n13, n14, n15;
   register double X, Y;
   double xy[15];
   mwIndex x0, y0, idx00; 

   /* node (0,0) */
   x0 = floor(x);
   y0 = floor(y);

   /* return NaN if too close to the edge */
   if (x0 < 2 || x0 >= nc-1 || y0 < 2 || y0 >= nr-1) {
      return NaN;
   }

   /* linear index of node (0,0) in map */
   idx00 = (x0-1)*nr + y0 -1;   
   z0  = pM[idx00];

   /* evalute the polynomial at (x,y) */
   X = x - x0;
   Y = y - y0;
   if ( (X == 0.0) && (Y == 0.0) )
      return z0;

   /* re-calculate coefficients only if necessary */
   if (!int_init || x0 != x0_last || y0 != y0_last) {
   
      z1  = pM[idx00 + nr];       /* node 1 at ( 1, 0) */
      z2  = pM[idx00 + nr + 1];   /* node 2 at ( 1, 1) */
      z3  = pM[idx00      + 1];   /* node 3 at ( 0, 1) */

      /* the 2nd ring nodes are needed for derivatives */
      n4  = pM[idx00 + tnr];     /* node 4 at  ( 2, 0) */
      n5  = pM[idx00 + tnr + 1]; /* node 5 at  ( 2, 1) */
      n6  = pM[idx00 + tnr + 2]; /* node 6 at  ( 2, 2) */
      n7  = pM[idx00 +  nr + 2]; /* node 7 at  ( 1, 2) */
      n8  = pM[idx00       + 2]; /* node 8 at  ( 0, 2) */
      n9  = pM[idx00 -  nr + 2]; /* node 9 at  (-1, 2) */
      n10 = pM[idx00 -  nr + 1]; /* node 10 at (-1, 1) */
      n11 = pM[idx00 -  nr];     /* node 11 at (-1, 0) */
      n12 = pM[idx00 -  nr - 1]; /* node 12 at (-1,-1) */
      n13 = pM[idx00       - 1]; /* node 13 at ( 0,-1) */
      n14 = pM[idx00 +  nr - 1]; /* node 14 at ( 1,-1) */
      n15 = pM[idx00 + tnr - 1]; /* node 15 at ( 2,-1) */

      /* the following "nodes" are derivatives at each point
       * in the x- or y-direction, calculated with central 
       * difference formulae 
       */
      z4  =  0.5 * (z1 - n11);             /* Dx@0 */
      z5  =  0.5 * (z3 - n13);             /* Dy@0 */   
      z6  =  0.5 * (n4 - z0);              /* Dx@1 */
      z7  =  0.5 * (z2 - n14);             /* Dy@1 */
      z8  =  0.5 * (n5 - z3);              /* Dx@2 */
      z9  =  0.5 * (n7 - z1);              /* Dy@2 */
      z10 =  0.5 * (z2 - n10);             /* Dx@3 */
      z11 =  0.5 * (n8 - z0);              /* Dy@3 */

      /*
       * Cross derivatives
       */ 
      z12 = 0.25 * (z2 - n10 - n14 + n12); /* Dxy@0 */
      z13 = 0.25 * (n5 - z3  - n15 + n13); /* Dxy@1 */
      z14 = 0.25 * (n6 - n8  - n4  + z0);  /* Dxy@2 */
      z15 = 0.25 * (n7 - n9  - z1  + n11); /* Dxy@3 */

      /* coefficients of the interpolation polynomial */
      ac[0]  = z4;
      ac[1]  = z5;
      ac[2]  = z12;
      ac[3]  = 3.0*(z1-z0)  -2.0*z4  -z6;
      ac[4]  = 3.0*(z3-z0)  -2.0*z5  -z11;
      ac[5]  = 3.0*(z10-z4) -2.0*z12 -z15;
      ac[6]  = 3.0*(z7-z5)  -2.0*z12 -z13;
      ac[7]  = 9.0*(z0-z1+z2-z3) +6.0*(z4+z5-z7-z10) +
	            3.0*(z6-z8-z9+z11) +4.0*z12 +2.0*(z13+z15) +z14;
      ac[8]  = 2.0*(z0-z1)  +z4  +z6;
      ac[9]  = 2.0*(z0-z3)  +z5  +z11;
      ac[10] = 2.0*(z4-z10) +z12 +z15;
      ac[11] = 2.0*(z5-z7)  +z12 +z13;
      ac[12] = 6.0*(z1-z0-z2+z3) +4.0*(z10-z4) +3.0*(z7-z5+z9-z11) +
               2.0*(z8-z6-z12-z15) -z13 -z14;
      ac[13] = 6.0*(z1-z0-z2+z3) +4.0*(z7-z5) +3.0*(z8-z4-z6+z10) +
               2.0*(z9-z11-z12-z13) -z14 -z15;
      ac[14] = 4.0*(z0-z1+z2-z3) +2.0*(z4+z5+z6-z7-z8-z9-z10+z11) +
               z12 +z13 +z14 +z15;

      int_init = 1;
   }

   /* store node 0 location */
   x0_last = x0;
   y0_last = y0;

   /* evaluate polynomial */ 
   xy[0]  = X;
   xy[1]  = Y;
   xy[2]  = X*Y;          /* x*y */
   xy[3]  = X*X;          /* x^2 */
   xy[4]  = Y*Y;          /* y^2 */
   xy[5]  = X*xy[4];      /* x*y^2 */
   xy[6]  = Y*xy[3];      /* y*x^2 */
   xy[7]  = xy[3]*xy[4];  /* x^2 * y^2 */
   xy[8]  = X*xy[3];      /* x^3 */ 
   xy[9]  = Y*xy[4];      /* y^3 */
   xy[10] = X*xy[9];      /* x*y^3 */
   xy[11] = Y*xy[8];      /* y*x^3 */
   xy[12] = X*xy[10];     /* x^2 * y^3 */
   xy[13] = Y*xy[11];     /* y^2 * x^3 */
   xy[14] = xy[8]*xy[9];  /* x^3 * y^3 */

   return z0 + dot_prod(15,ac,xy); 
}


/*-----------------------------------------------------------------*/

static double 
interp_pchip(double * RESTRICT pM, double x, double y) 
{
   double z0, z1, z2, z3, z4, z5, z6, z7, z8, z9, z10, z11, z12, z13, z14, z15;
   double b0, b1, b2, b3;
   register double T;
   double A[3];
   double t[3];
   mwIndex x0, y0, idx00; 

   /* node (0,0) */
   x0 = floor(x);
   y0 = floor(y);

   /* return NaN if too close to the edge */      
   if (x0 < 2 || x0 >= nc-1 || y0 < 2 || y0 >= nr-1) {
      return NaN;
   }

   /* matrix values at nodes */
   idx00 = (x0-1)*nr +y0 -1;   /* linear index of node (0,0) in map */
   z0  = pM[idx00];            /* node 0 at  ( 0, 0) */

   /* evalute the polynomial at (x,y) */
   if ( (x - x0 == 0.0) && (y - y0 == 0.0) )
      return z0;

   z1  = pM[idx00 + nr];       /* node 1 at  ( 1, 0) */
   z2  = pM[idx00 + nr + 1];   /* node 2 at  ( 1, 1) */
   z3  = pM[idx00      + 1];   /* node 3 at  ( 0, 1) */
   z4  = pM[idx00 + tnr];      /* node 4 at  ( 2, 0) */
   z5  = pM[idx00 + tnr + 1];  /* node 5 at  ( 2, 1) */
   z6  = pM[idx00 + tnr + 2];  /* node 6 at  ( 2, 2) */
   z7  = pM[idx00 +  nr + 2];  /* node 7 at  ( 1, 2) */
   z8  = pM[idx00       + 2];  /* node 8 at  ( 0, 2) */
   z9  = pM[idx00 -  nr + 2];  /* node 9 at  (-1, 2) */
   z10 = pM[idx00 -  nr + 1];  /* node 10 at (-1, 1) */
   z11 = pM[idx00 -  nr];      /* node 11 at (-1, 0) */
   z12 = pM[idx00 -  nr - 1];  /* node 12 at (-1,-1) */
   z13 = pM[idx00       - 1];  /* node 13 at ( 0,-1) */
   z14 = pM[idx00 +  nr - 1];  /* node 14 at ( 1,-1) */
   z15 = pM[idx00 + tnr - 1];  /* node 15 at ( 2,-1) */

   /* interpolation in x direction */
   T = x - x0;
   t[0] = T;
   t[1] = T*T;
   t[2] = T*t[1];
   
   /* b0 from z12, z13, z14, z15 */
   A[0] = z14 - z12;
   A[1] = 2.0*z12 - 5.0*z13 + 4.0*z14 - z15;
   A[2] = 3.0*(z13-z14) - z12 + z15;
   b0 = z13 + 0.5 * dot_prod(3, A, t);
   
   /* b1 from z11, z0, z1, z4 */
   A[0] = z1 - z11;
   A[1] = 2.0*z11 - 5.0*z0 + 4.0*z1 - z4;
   A[2] = 3.0*(z0-z1) - z11 + z4;
   b1 = z0 + 0.5 * dot_prod(3, A, t);
   
   /* b2 from z10, z3, z2, z5 */
   A[0] = z2 - z10;
   A[1] = 2.0*z10 - 5.0*z3 + 4.0*z2 - z5;
   A[2] = 3.0*(z3-z2) - z10 + z5;
   b2 = z3 + 0.5 * dot_prod(3, A, t);
   
   /* b3 from z9, z8, z7, z6 */
   A[0] = z7 - z9;
   A[1] = 2.0*z9 - 5.0*z8 + 4.0*z7 - z6;
   A[2] = 3.0*(z8-z7) - z9 + z6;
   b3 = z8 + 0.5 * dot_prod(3, A, t);
   
   /* interpolation in y direction */
   T = y - y0;
   t[0] = T;
   t[1] = T*T;
   t[2] = T*t[1];
   
   A[0] = b2 - b0;
   A[1] = 2.0*b0 - 5.0*b1 + 4.0*b2 - b3;
   A[2] = 3.0*(b1-b2) - b0 + b3;

   return b1 + 0.5 * dot_prod(3, A, t);
}


/*-----------------------------------------------------------------*/

static double 
interp_lagrange(double * RESTRICT pM, double x, double y)
{
   register double iden, X, Y;
   double z[4];
   double c[4];
   double d[4];
   double w[4] = {-1.0, 3.0, -3.0, 1.0};
   mwIndex x0, y0, idx00; 

   /* node (0,0) */
   x0 = floor(x);
   y0 = floor(y);

   /* return NaN if too close to the edge */      
   if (x0 < 2 || x0 >= nc-1 || y0 < 2 || y0 >= nr-1) {
      return NaN;
   }

   /* linear index of node (0,0) in map */
   idx00 = (x0-1)*nr + y0 -1;  

   /* evalute the polynomial at (x,y) */
   X = x - x0;
   Y = y - y0;

   /* x-interpolation */
   if ( X == 0.0 ) { /* point is at x0 */
      c[0] = pM[idx00  - 1];  /* node 13 at ( 0,-1) */
      c[1] = pM[idx00];       /* node 0 at  ( 0, 0) */
      c[2] = pM[idx00  + 1];  /* node 3 at  ( 0, 1) */
      c[3] = pM[idx00  + 2];  /* node 8 at  ( 0, 2) */
   }
   else {
      d[0] = 1.0/(X + 1);
      d[1] = 1.0/X;
      d[2] = 1.0/(X - 1);
      d[3] = 1.0/(X - 2);
      iden = 1.0 / dot_prod(4,d,w);

      z[0] = pM[idx00 -  nr - 1];  /* node 12 at (-1,-1) */
      z[1] = pM[idx00       - 1];  /* node 13 at ( 0,-1) */
      z[2] = pM[idx00 +  nr - 1];  /* node 14 at ( 1,-1) */
      z[3] = pM[idx00 + tnr - 1];  /* node 15 at ( 2,-1) */
      c[0] = iden * dot_prodw(4,z,d,w);

      z[0] = pM[idx00 - nr];       /* node 11 at (-1, 0) */
      z[1] = pM[idx00];            /* node 0 at  ( 0, 0) */
      z[2] = pM[idx00 + nr];       /* node 1 at  ( 1, 0) */
      z[3] = pM[idx00 + tnr];      /* node 4 at  ( 2, 0) */
      c[1] = iden * dot_prodw(4,z,d,w);

      z[0] = pM[idx00 -  nr + 1];  /* node 10 at (-1, 1) */
      z[1] = pM[idx00       + 1];  /* node 3 at  ( 0, 1) */
      z[2] = pM[idx00 + nr  + 1];  /* node 2 at  ( 1, 1) */
      z[3] = pM[idx00 + tnr + 1];  /* node 5 at  ( 2, 1) */
      c[2] = iden * dot_prodw(4,z,d,w);

      z[0] = pM[idx00 -  nr + 2];  /* node 9 at  (-1, 2) */
      z[1] = pM[idx00       + 2];  /* node 8 at  ( 0, 2) */
      z[2] = pM[idx00 +  nr + 2];  /* node 7 at  ( 1, 2) */
      z[3] = pM[idx00 + tnr + 2];  /* node 6 at  ( 2, 2) */
      c[3] = iden * dot_prodw(4,z,d,w);
   }
   
   /* y-interpolation */
   if ( Y == 0.0 ) {
      return c[1];
   }
   else {
      d[0] = 1.0/(Y + 1);
      d[1] = 1.0/Y;
      d[2] = 1.0/(Y - 1);
      d[3] = 1.0/(Y - 2);
      return dot_prodw(4,d,c,w) / dot_prod(4,d,w);
   }
}


/*-----------------------------------------------------------------*/

static INLINE double 
dot_prod(int n, double * RESTRICT a, double * RESTRICT b) 
{
   register int k;
   register double sum = 0.0;

   for (k=0; k<n; k++) 
      sum += a[k]*b[k];

   return sum;
}


/*-----------------------------------------------------------------*/

static INLINE double 
dot_prodw(int n, double * RESTRICT a, 
	         double * RESTRICT b, double * RESTRICT w) 
{
   register int k;
   register double sum = 0.0;

   for (k=0; k<n; k++) 
      sum += a[k]*b[k]*w[k];

   return sum;
}

/*-- 终 -------------------------------------------------------*/
