 /*
 * =============================================================
 * mexBasinModel.c - C MEX file to mark drainage basins
 *
 * usage: [Basin, A] = mexBasinModel(D,i,j,bvec), where D is a K x J matrix of
 *        D8 drainage directions, i,j are vectors of indices of basin outlets, and
 *        bvec is a vector of codes for the [left, right, upper, lower]
 *        boundary conditions (0=fixed, 1=mirror, 2=periodic)
 *
 *        Returns Basin, a K x J matrix of integers marking basins
 *        that drain to the points in [i,j], and A, D8 TCA in cells.
 *
 * This is a MEX-file for MATLAB.
 * Copyright (c) 2008-2011 Taylor Perron
 * =============================================================
 */

#include "mex.h"
#include "matrix.h"
#include <math.h>
#include <time.h>
#include <stdlib.h>

#define FIXED    0
#define MIRROR   1
#define PERIODIC 2

#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )


void GetArea(const int i, const int j, const int K, const int J, const int basinnum, double A[], double Basin[], double D[], double bdy[])
{

int k, p, q, c=K*j+i, bl, br, bu, bd;
int r[] = {5,6,7,8,1,2,3,4};
int di[] = {0,-1,-1,-1,0,1,1,1};
int dj[] = {1,1,0,-1,-1,-1,0,1};
int skip[] = {0,0,0,0,0,0,0,0};
int w[] = {1,1,1,1,1,1,1,1};

if (A[c]>0) {return;} // Bail out if we already know the area

A[c] = 1; /* each cell drains at least itself */

Basin[c] = basinnum; /* mark this cell as part of the present basin */

// left
if (bdy[0] == 0) {
    bl=FIXED;
} else if (bdy[0] == 1) {
    bl=MIRROR;
} else if (bdy[0] == 2) {
    bl=PERIODIC;
}

// right
if (bdy[1] == 0) {
    br=FIXED;
} else if (bdy[1] == 1) {
    br=MIRROR;
} else if (bdy[1] == 2) {
    br=PERIODIC;
}

// upper
if (bdy[2] == 0) {
    bu=FIXED;
} else if (bdy[2] == 1) {
    bu=MIRROR;
} else if (bdy[2] == 2) {
    bu=PERIODIC;
}

// lower
if (bdy[3] == 0) {
    bd=FIXED;
} else if (bdy[3] == 1) {
    bd=MIRROR;
} else if (bdy[3] == 2) {
    bd=PERIODIC;
}


if (i==0) {
    
    switch (bu) {
    case FIXED: // fixed upper
        skip[1]=1;
        skip[2]=1;
        skip[3]=1;
        break;
    case MIRROR: // mirror upper
        skip[1]=1;
        skip[2]=1;
        skip[3]=1;
        w[5] = 2;
        w[6] = 2;
        w[7] = 2;
        break;
    case PERIODIC: // periodic upper
        di[1]=K-1;
        di[2]=K-1;
        di[3]=K-1;
        break;
    }

} else if (i==K-1) {

    switch (bd) {
    case FIXED: // fixed lower
        skip[5]=1;
        skip[6]=1;
        skip[7]=1;
        break;
    case MIRROR: // mirror lower
        skip[5]=1;
        skip[6]=1;
        skip[7]=1;
        w[1] = 2;
        w[2] = 2;
        w[3] = 2;
        break;
    case PERIODIC: // periodic lower
        di[5]=1-K;
        di[6]=1-K;
        di[7]=1-K;
        break;
    }


}

if (j==0) {
    
    switch (bl) {
    case FIXED: // fixed left
        skip[3]=1;
        skip[4]=1;
        skip[5]=1;
        break;
    case MIRROR: // mirror left
        skip[3]=1;
        skip[4]=1;
        skip[5]=1;
        w[7] = 2;
        w[0] = 2;
        w[1] = 2;
        break;
    case PERIODIC: // periodic left
        dj[3]=J-1;
        dj[4]=J-1;
        dj[5]=J-1;
        break;
    }
    
} else if (j==J-1) {

    switch (br) {
    case FIXED: // fixed right
        skip[7]=1;
        skip[0]=1;
        skip[1]=1;
        break;
    case MIRROR: // mirror right
        skip[7]=1;
        skip[0]=1;
        skip[1]=1;
        w[3] = 2;
        w[4] = 2;
        w[5] = 2;
        break;
    case PERIODIC: // periodic right
        dj[7]=1-J;
        dj[0]=1-J;
        dj[1]=1-J;
        break;
    }

}


for (k=0; k<8; k++) {  /* loop through drainage directions */
    
    if (!skip[k]) {
        
        p = i + di[k];
        q = j + dj[k];
                
        if ((int)D[K*q+p] == r[k]) {  /* if the current drainage direction has a weight */
            
            GetArea(p, q, K, J, basinnum, A, Basin, D, bdy);  /* recursive call to get drainage area for neighbor k */
            A[c] += w[k] * A[K*q+p]; /* increment A(i,j) by the area of k */
            /* Purpose of the w prefactor: A point on a mirrored boundary receiving flow from a point off the boundary should receive twice the area */               
        }
    }
    
} /* end of drainage direction loop */

} /* end GetArea subfunction */



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

  double *M, *D, *F, *R, *O, *W, *B, *A, *Basin, *minima, *a, *i0, *j0, *bdy;
  int i, j, k, K, J, *numU, *numMin, numOut, outletidx[2], ndims=3, dims[]={0,0,8};

  // Get pointers to inputs
  D = (double *)mxGetPr(prhs[0]); // D8 directions
  i0 = (double *)mxGetPr(prhs[1]); // outlet i's
  j0 = (double *)mxGetPr(prhs[2]); // outlet j's
  bdy = (double *)mxGetPr(prhs[3]); // boundary condition vector
  
  // Get number of outlets
  numOut=max(mxGetM(prhs[1]),mxGetN(prhs[1])); // we do it this way b/c we don't know if they're row or col vectors
  
//   // convert outlet indices from 1-based to zero-based
//   *i0 = *i0-1;
//   *j0 = *j0-1;
  
  // Get dimensions of input matrix of directions
  K=dims[0]=mxGetM(prhs[0]);
  J=dims[1]=mxGetN(prhs[0]);


  // Create array for the return arguments  
  Basin = (double *)mxGetPr(plhs[0]= mxCreateDoubleMatrix(K, J, mxREAL)); // integers labeling basins that drain to [i0,j0], zero otherwise
  A = (double *)mxGetPr(plhs[1]= mxCreateDoubleMatrix(K, J, mxREAL)); // TCA for cells that contribute flow to i0,j0, zero otherwise



  ////////////////////////////////////////
  // Now run  GetArea routine to get TCA
  ////////////////////////////////////////

  // Get area only for cells that contribute some flow to the specified cells
  // Note that we convert from 1-based to zero-based indices on the fly
  for (k=0; k<numOut; k++) {
    GetArea((int)i0[k]-1,(int)j0[k]-1,K,J,k+1,A,Basin,D,bdy);
  }

} // End mexFunction
