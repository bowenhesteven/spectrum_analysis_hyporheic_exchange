/*
* ========================================================================
* mexD8.c - C MEX file to compute upslope contributing area using the D8
*           (steepest descent) algorithm of O'Callaghan & Mark (1984)
*
*     usage: [A, B, D] = mexD8(M,dy/dx,flood,ndv)
*
*     Input arguments:
*       M      a K x J matrix of elevations
*       dy/dx  the ratio of grid spacings in the y and x directions. e.g.,
*              if dy=4m and dx=2m, dy/dx=2. (optional: default = 1)
*       flood  if 1, routes flow through local minima in M (potentially
*              time-consuming). if zero, skips this step (faster).
*              (optional: default = 0)
*       ndv    matrix value that indicates no data. this will be used to
*              identify the boundaries of your elevation matrix. If all
*              the elements in your matrix have valid elevations, just set
*              ndv to some value that does not appear in your elevations
*              (e.g., -1 or 9999). (optional: default = 0)
*
*     Return arguments:
*       A      a matrix of total contributing areas (in units of CELLS)
*       B      a matrix of cells that receive drainage from boundaries
*       D      a matrix of D8 drainage directions (1-8, CCW from east, zero
*              at sinks and boundaries) 
*
* This is a MEX-file for MATLAB.
*
* Copyright (C) 2004-2010 Taylor Perron <perron@mit.edu>
* 
* This program is free software: you can redistribute it and/or modify it 
* under the terms of the GNU General Public License as published by the 
* Free Software Foundation. You should have received a copy of the GNU 
* General Public License along with this program.  If not, see 
* http://www.gnu.org/licenses.
* ========================================================================
*/

#include "mex.h"
#include "matrix.h"
#include <math.h>
  
void GetBound(const int K, const int J, double M[], double Bdy[], double B[], double *ndv)
{

int i, j, Kj, Kju, Kjd;
double n;

n = *ndv;

for (j=0; j<J; j++) {
    Kj=K*j;
    Kju=K*(j+1);
    Kjd=K*(j-1);
    for (i=0; i<K; i++) {
        if (M[Kj+i] != n) { // if M(i,j) has data
            if (i==0 || i==K-1 || j==0 || j==J-1) { // if on matrix boundary
                Bdy[Kj+i] = 2;
                B[Kj+i] = 1; // mark as boundary-influenced
            } else if (M[Kju+i] == n || M[Kju+i-1] == n || M[Kj+i-1] == n || M[Kjd+i-1] == n || M[Kjd+i] == n || M[Kjd+i+1] == n || M[Kj+i+1] == n || M[Kju+i+1] == n) { // if any of its neighbors has nodata
                Bdy[Kj+i] = 2;
                B[Kj+i] = 1; // mark as boundary-influenced
            } else { // otherwise it's an interior elevation cell
                Bdy[Kj+i] = 1;
            }
        } // if M(i,j)==n, leave Bdy=0
    }
}
    
} // end GetBound()


void D8Dir(const int K, const int J, double M[], double Bdy[], double D[], double minima[], double *a, int *numMin)
{

int i, j, k, Kj, Kjup, Kjdown, theD;
double invdx, invdy, invdiag, theS, maxS;
double e0;

// For boundaries and corners, restrict facets to the following:
int numfacets, *facetidx;
int ULfacets[]={6,7,0};
int URfacets[]={4,5,6};
int LLfacets[]={0,1,2};
int LRfacets[]={2,3,4};
int Lfacets[]={6,7,0,1,2};
int Rfacets[]={2,3,4,5,6};
int Tfacets[]={4,5,6,7,0};
int Bfacets[]={0,1,2,3,4};
int Allfacets[]={0,1,2,3,4,5,6,7};

invdx=1;
invdy=1/(*a);
invdiag=1/sqrt(1+(*a)*(*a));

  for (j=0; j<J; j++) { // Loop through columns

      Kj=K*j;
      Kjup=K*(j+1);
      Kjdown=K*(j-1);

      for (i=0; i<K; i++) { // Loop through rows
          
          if (Bdy[Kj+i] > 0) { // if it is nodata, skip it

              e0=M[Kj+i]; // Elevation at (i,j)
              maxS = 0; // Set max slope to zero
              theD = 0;

              if (i>0 && i<(K-1) && j>0 && j<(J-1)) { // Interior cell
                  numfacets=8;
                  facetidx=Allfacets;
              } else if (i==0 && j==0) { // UL corner
                  numfacets=3;
                  facetidx=ULfacets;
              } else if (i==0 && j==J-1) { // UR corner
                  numfacets=3;
                  facetidx=URfacets;
              } else if (i==K-1 && j==0) { // LL corner
                  numfacets=3;
                  facetidx=LLfacets;
              } else if (i==K-1 && j==J-1) { // LR corner
                  numfacets=3;
                  facetidx=LRfacets;
              } else if (j==0) { // Left boundary
                  numfacets=5;
                  facetidx=Lfacets;
              } else if (j==J-1) { // Right boundary
                  numfacets=5;
                  facetidx=Rfacets;
              } else if (i==0) { // Top boundary
                  numfacets=5;
                  facetidx=Tfacets;
              } else if (i==K-1) { // Bottom boundary
                  numfacets=5;
                  facetidx=Bfacets;
              }

              for (k=0; k<numfacets; k++) { // Loop through facets, ignoring any that have nodata
                  switch (facetidx[k]) {
                      case 0:
                        theS=e0-M[Kjup+i];
                        if (Bdy[Kjup+i] == 0) {
                            theS = -1;
                        }
                        break;
                      case 1:
                        theS=(e0-M[Kjup+(i-1)])*invdiag;
                        if (Bdy[Kjup+(i-1)] == 0) {
                            theS = -1;
                        }
                        break;
                      case 2:
                        theS=(e0-M[Kj+(i-1)])*invdy;
                        if (Bdy[Kj+(i-1)] == 0) {
                            theS = -1;
                        }
                        break;
                      case 3:
                        theS=(e0-M[Kjdown+(i-1)])*invdiag;
                        if (Bdy[Kjdown+(i-1)] == 0) {
                            theS = -1;
                        }
                        break;
                      case 4:
                        theS=e0-M[Kjdown+i];
                        if (Bdy[Kjdown+i] == 0) {
                            theS = -1;
                        }
                        break;
                      case 5:
                        theS=(e0-M[Kjdown+(i+1)])*invdiag;
                        if (Bdy[Kjdown+(i+1)] == 0) {
                            theS = -1;
                        }
                        break;
                      case 6:
                        theS=(e0-M[Kj+(i+1)])*invdy;
                        if (Bdy[Kj+(i+1)] == 0) {
                            theS = -1;
                        }
                        break;
                      case 7:
                        theS=(e0-M[Kjup+(i+1)])*invdiag;
                        if (Bdy[Kjup+(i+1)] == 0) {
                            theS = -1;
                        }
                        break;
                  }

                  // If this is a downslope facet AND the steepest yet, record its slope and direction
                  if (theS>maxS) {
                      maxS=theS;
                      theD=facetidx[k]+1;
                  }

              } // end looping through facets

              // If (i,j) is a local minimum (maxS remains zero), it retains a D8 direction of zero and is flagged in minima
              if (maxS==0) {
                  minima[Kj+i]=1;
                  (*numMin)++;
              }
              D[Kj+i]=theD; // assign D8 direction
          }
      } 
  } 
} // end D8Dir()


void D8Resolve(const int i, const int j, const int K, const int J, double D[], double F[], double R[], int *numU)
{

int k, p, q, r;

if (R[K*j+i]<1) { // Only proceed if this cell is not yet resolved

    R[K*j+i]=1; // Now it is resolved
    (*numU)--; // One less cell to resolve

        for (k=1; k<9; k++) {

            switch (k) {
              case 1:
                p=i; q=j+1; r=5;
                break;
              case 2:
                p=i-1; q=j+1; r=6;
                break;
              case 3:
                p=i-1; q=j; r=7;
                break;
              case 4:
                p=i-1; q=j-1; r=8;
                break;
              case 5:
                p=i; q=j-1; r=1;
                break;
              case 6:
                p=i+1; q=j-1; r=2;
                break;
              case 7:
                p=i+1; q=j; r=3;
                break;
              case 8:
                p=i+1; q=j+1; r=4;
                break;
            }



    	if (p>=0 && p<K && q>=0 && q<J) { // avoid going off boundaries

            if (D[K*q+p] == r) { // if neighbor k drains to current cell
                D8Resolve(p,q,K,J,D,F,R,numU);
            }

        } // if statement that avoids going off boundaries
    } // drainage direction loop
}

} // End D8Resolve()


void FindOutlet(const int K, const int J, double M[], double R[], double O[], int outletidx[])
{

int i, j, jdown, jup, Kj, Kjdown, Kjup, nsum;
double outletelev=1E20; // A larger elevation than we are likely to encounter

for (j=1; j<(J-1); j++) { // Loop through columns, excluding x boundaries

  jdown = j-1; jup = j+1;

  Kj=K*j;
  Kjup=K*jup;
  Kjdown=K*jdown;

  for (i=1; i<(K-1); i++) { // Loop through rows, excluding y boundaries

      if (R[K*j+i] && !O[K*j+i]) { // If it's not resolved, or has already been used as an outlet, don't bother with it
          nsum = R[Kjup+i] + R[Kjup+(i-1)] + R[Kj+(i-1)] + R[Kjdown+(i-1)] + R[Kjdown+i] + R[Kjdown+i+1] + R[Kj+i+1] + R[Kjup+i+1];
          if (nsum<8) {
              if (M[K*j+i] < outletelev) {
                  outletelev = M[K*j+i];
                  outletidx[0] = i;
                  outletidx[1] = j;
              }
          }
      }

  }

  i=0; // Upper y boundary
      if (R[K*j+i] && !O[K*j+i]) {
          // Only go through neighbors 5-1 for upper boundary
          nsum = R[Kjup+i] + R[Kjdown+i] + R[Kjdown+i+1] + R[Kj+i+1] + R[Kjup+i+1];
          if (nsum<5) {
              if (M[K*j+i] < outletelev) {
                  outletelev = M[K*j+i];
                  outletidx[0] = i;
                  outletidx[1] = j;
              }
          }
      }

  i=K-1; // Lower y boundary
      if (R[K*j+i] && !O[K*j+i]) {
          // Only go through neighbors 1-5 for lower boundary
          nsum = R[Kjup+i] + R[Kjup+(i-1)] + R[Kj+(i-1)] + R[Kjdown+(i-1)] + R[Kjdown+i];
          if (nsum<5) {
              if (M[K*j+i] < outletelev) {
                  outletelev = M[K*j+i];
                  outletidx[0] = i;
                  outletidx[1] = j;
              }
          }
      }

}

j=0; // Left x boundary

  jup = j+1;

  Kj=K*j;
  Kjup=K*jup;

  for (i=1; i<(K-1); i++) { // Loop through rows, excluding y boundaries

      if (R[K*j+i] && !O[K*j+i]) { // If it's not resolved, or has already been used as an outlet, don't bother with it
          nsum = R[Kjup+i] + R[Kjup+(i-1)] + R[Kj+(i-1)] + R[Kj+i+1] + R[Kjup+i+1];
          if (nsum<5) {
              if (M[K*j+i] < outletelev) {
                  outletelev = M[K*j+i];
                  outletidx[0] = i;
                  outletidx[1] = j;
              }
          }
      }

  }

  i=0; // Upper left corner
      if (R[K*j+i] && !O[K*j+i]) {
          // Only go through neighbors 7-1
          nsum = R[Kjup+i] + R[Kj+i+1] + R[Kjup+i+1];
          if (nsum<3) {
              if (M[K*j+i] < outletelev) {
                  outletelev = M[K*j+i];
                  outletidx[0] = i;
                  outletidx[1] = j;
              }
          }
      }

  i=K-1; // Lower left corner
      if (R[K*j+i] && !O[K*j+i]) {
          // Only go through neighbors 1-3
          nsum = R[Kjup+i] + R[Kjup+(i-1)] + R[Kj+(i-1)];
          if (nsum<3) {
              if (M[K*j+i] < outletelev) {
                  outletelev = M[K*j+i];
                  outletidx[0] = i;
                  outletidx[1] = j;
              }
          }
      }

j=J-1; // Right x boundary

  jdown = j-1;

  Kj=K*j;
  Kjdown=K*jdown;

  for (i=1; i<(K-1); i++) { // Loop through rows, excluding y boundaries

      if (R[K*j+i] && !O[K*j+i]) { // If it's not resolved, or has already been used as an outlet, don't bother with it
          nsum = R[Kj+(i-1)] + R[Kjdown+(i-1)] + R[Kjdown+i] + R[Kjdown+i+1] + R[Kj+i+1];
          if (nsum<5) {
              if (M[K*j+i] < outletelev) {
                  outletelev = M[K*j+i];
                  outletidx[0] = i;
                  outletidx[1] = j;
              }
          }
      }


  }

  i=0; // Upper right corner
      if (R[K*j+i] && !O[K*j+i]) {
          // Only go through neighbors 5-7
          nsum = R[Kjdown+i] + R[Kjdown+i+1] + R[Kj+i+1];
          if (nsum<3) {
              if (M[K*j+i] < outletelev) {
                  outletelev = M[K*j+i];
                  outletidx[0] = i;
                  outletidx[1] = j;
              }
          }
      }

  i=K-1; // Lower right corner
      if (R[K*j+i] && !O[K*j+i]) {
          // Only go through neighbors 3-5 for lower boundary
          nsum = R[Kj+(i-1)] + R[Kjdown+(i-1)] + R[Kjdown+i];
          if (nsum<3) {
              if (M[K*j+i] < outletelev) {
                  outletelev = M[K*j+i];
                  outletidx[0] = i;
                  outletidx[1] = j;
              }
          }
      }


} // End FindOutlet()


void Flood(const int K, const int J, int i, int j, double M[], double Bdy[], double D[], double R[], double F[], double waterlevel, int *numU)
{

    int k, p, q, r;

    for (k=0; k<8; k++) {

        switch (k+1) {
          case 1:
            p=i; q=j+1; r=5; // E
            break;
          case 2:
            p=i-1; q=j+1; r=6; // NE
            break;
          case 3:
            p=i-1; q=j; r=7; // N
            break;
          case 4:
            p=i-1; q=j-1; r=8; // NW
            break;
          case 5:
            p=i; q=j-1; r=1; // W
            break;
          case 6:
            p=i+1; q=j-1; r=2; // SW
            break;
          case 7:
            p=i+1; q=j; r=3; // S
            break;
          case 8:
            p=i+1; q=j+1; r=4; // SE
            break;
        }


        if (p>=0 && p<K && q>=0 && q<J && Bdy[K*q+p]>0) { // avoid going off boundaries or into nodata regions

            if (!R[K*q+p]) {
                // if neighbor k is at or below the water level
                if (M[K*q+p] <= waterlevel) { // neighbor k has equal or lower elevation

                    F[K*q+p] = 1; // we flooded it
                    R[K*q+p] = 1; // B/C current cell is resolved (drains to a boundary), neighbor k is now resolved
                    (*numU)--;  // one less cell to worry about
                    D[K*q+p] = r; // make neighbor k drain toward current cell
                    Flood(K,J,p,q,M,Bdy,D,R,F,waterlevel,numU); // Call flood recursively for neighbor k

                } else if (D[K*q+p]==r) { // it's above the water level and it drains to the current cell

                    D8Resolve(p,q,K,J,D,F,R,numU);

                } // otherwise, it has a higher elevation but drains to another as yet unresolved cell; leave it alone
            } // if statement that checks if (p,q) is already resolved
        } // if statement that avoids going off y boundaries
    } // drainage direction loop


} // end Flood()



void GetD8Area(const int i, const int j, const int K, const int J, double A[], double D[], double B[], double Bdy[])
{

int k, p, q, r, Kqplusp, c=K*j+i;

if (A[c]>0) {return;} // Bail out if we already know the area

A[c] = 1; // each cell drains at least itself

for (k=0; k<8; k++) {  // loop through drainage directions

    switch (k) {
      case 0:
        p=i; q=j+1; r=5;
        break;
      case 1:
        p=i-1; q=j+1; r=6;
        break;
      case 2:
        p=i-1; q=j; r=7;
        break;
      case 3:
        p=i-1; q=j-1; r=8;
        break;
      case 4:
        p=i; q=j-1; r=1;
        break;
      case 5:
        p=i+1; q=j-1; r=2;
        break;
      case 6:
        p=i+1; q=j; r=3;
        break;
      case 7:
        p=i+1; q=j+1; r=4;
        break;
    }

    if (p>=0 && p<K && q>=0 && q<J && Bdy[K*q+p]>0) {  // avoid going off boundaries or into nodata regions

        Kqplusp = K*q+p;

        if (D[Kqplusp] == r) {  // if neighbor k drains to current cell
            GetD8Area(p,q,K,J,A,D,B,Bdy);  // recursive call to get drainage area for neighbor k
            A[c] += A[Kqplusp]; // increment A(i,j) by the area of k

            // Now we know whether any of the boundaries drain to neighbor p,q (or if it is a boundary itself).
            // If so, the boundaries drain to i,j as well, so set B(i,j)=1
            if (B[Kqplusp]==1) {
                B[c]=1;
            }

        }
    } // end of if statement that avoids going off boundaries 
} // end of drainage direction loop 

} // end GetD8Area()



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

  mxArray *Rptr, *Optr;
  double *M, *Bdy, *a, *fl, *ndv, *A, *B, *D, *F, *R, *O, *minima;
  int i, j, K, J, *numU, *numMin, outletidx[2], ndims=3, dims[]={0,0,8};

  // Argument checking
  if (nrhs < 1 || nrhs > 4) {
    mexErrMsgTxt("Wrong number of input arguments. Usage: [area boundary] = mexD8(elev,dy/dx,doflooding,nodataval)");
  } else if (nlhs > 3) {
    mexErrMsgTxt("Too many output arguments. Usage: [area boundary] = mexD8(elev,dy/dx,doflooding,nodataval)");
  }

  if (mxGetNumberOfDimensions(prhs[0]) != 2) {
    mexErrMsgTxt("Elevation array must be 2D.");
  }
  
  // Get pointers to inputs, assign default values if not supplied
  M = (double *)mxGetPr(prhs[0]); // elevations

  if (nrhs < 2) { 
      a = (double *)mxGetPr(mxCreateDoubleMatrix(1, 1, mxREAL));
      *a = 1; // default
  } else {
      a = (double *)mxGetPr(prhs[1]); // dy/dx
  }
  
  if (nrhs < 3) { 
      fl = (double *)mxGetPr(mxCreateDoubleMatrix(1, 1, mxREAL));
      *fl = 0; // default
  } else {
      fl = (double *)mxGetPr(prhs[2]); // flood local minima, or not
  }

  if (nrhs < 4) { 
      ndv = (double *)mxGetPr(mxCreateDoubleMatrix(1, 1, mxREAL));
      *ndv = 0; // default
  } else {
      ndv = (double *)mxGetPr(prhs[3]); // nodata value
  }


  // Get dimensions of input matrix of elevations
  K=dims[0]=mxGetM(prhs[0]);
  J=dims[1]=mxGetN(prhs[0]);

  // Create arrays for the return arguments
  A = (double *)mxGetPr(plhs[0]= mxCreateDoubleMatrix(K, J, mxREAL)); // matrix of TCAs
  B = (double *)mxGetPr(plhs[1]= mxCreateDoubleMatrix(K, J, mxREAL)); // indicates which cells receive drainage from one or more boundary cells
  D = (double *)mxGetPr(plhs[2]= mxCreateDoubleMatrix(K, J, mxREAL)); // matrix of D8 drainage directions

  // Create internally used arrays
  Bdy = (double *)mxGetPr(mxCreateDoubleMatrix(K, J, mxREAL)); // matrix of boundary flags
//   D = (double *)mxGetPr(mxCreateDoubleMatrix(K, J, mxREAL)); // matrix of D8 drainage directions
  minima = (double *)mxGetPr(mxCreateDoubleMatrix(K, J, mxREAL)); // 1 if an element is a local minimum, zero otherwise
  F = (double *)mxGetPr(mxCreateDoubleMatrix(K, J, mxREAL)); // matrix that indicates what was flooded
  numMin = mxCalloc(1, sizeof(int)); // number of local minima

  // Find boundary elements, defined as cells with valid elevations that
  // are on the edge of the grid or have at least one neighbor with nodata. 
  // Bdy == 1 for interior points, 2 on boundaries, zero for nodata
  GetBound(K,J,M,Bdy,B,ndv);

  // Calculate D8 drainage directions
  *numMin = 0;
  D8Dir(K,J,M,Bdy,D,minima,a,numMin); // Calculate D8 drainage directions

  // if there are local minima and the user opted to route flow through
  // minima, do the flooding routines  
  if (*numMin > 0 && *fl == 1) {
      
      numU = mxCalloc(1, sizeof(int)); // number of cells with unresolved drainage
      *numU = K*J; // initially, all cells are considered unresolved

      // Create internally used arrays
      R = (double *)mxGetPr(Rptr= mxCreateDoubleMatrix(K, J, mxREAL)); // indicates which cells have resolved drainage
      O = (double *)mxGetPr(Optr= mxCreateDoubleMatrix(K, J, mxREAL)); // indicates which cells have been used as outlets

      for (i=0; i<K; i++) { // Find cells that drain to one of the data boundaries
          for (j=0; j<J; j++) {
              if (Bdy[K*j+i] == 2) {
                  D8Resolve(i,j,K,J,D,F,R,numU);
              } else if (Bdy[K*j+i] == 0) {
                  (*numU)--; // remove the nodata cells from the count of unresolved cells
              }
          }
      }

      while ((*numU)>0) {
          FindOutlet(K,J,M,R,O,outletidx); 
          // start at outletidx and resolve "upslope" using Flood, decreasing numU and recording used outlets in O[] as you go
          O[K*outletidx[1]+outletidx[0]] = 1; // We're about to use this as an outlet
          Flood(K,J,outletidx[0],outletidx[1],M,Bdy,D,R,F,M[K*outletidx[1]+outletidx[0]],numU);
      }

      // free memory
      mxFree(numU);
      mxFree(numMin);
      mxDestroyArray(Rptr);
      mxDestroyArray(Optr);

  } // end of flooding routines

  
  // Get contributing area for each element
  for (i=0; i<K; i++) {
    for (j=0; j<J; j++) {
      GetD8Area(i,j,K,J,A,D,B,Bdy); // D8
    }
  }
  
} // end mexFunction()
