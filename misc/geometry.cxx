#include "../misc/constants.h"

extern double **ptr;
extern int **ntr,**nne;

/*
  Input parameters:
  nge - Geometry type
  rb  - Rough radius of domain
  ndv - Number of wall element divisions

  Output parameters:
  nwe - Number of elements
  nwn - Number of nodes
  nwc - Unused
  ptr - Coordinates
  ntr - Connectivity
  nne - Dual graph
*/
void geometry(int ndv, int& nwc, int nge, double rb) {
  // nen - Number of nodes per element
  int i,k,j,iz,nen,l,nm,ic;
  int n1[4],n2[4],i2[2],i3[2];
  double rad,radc;
  // Coordinates
  //   The first three entries are the nodes at triangle vertices
  //   The second three are midnodes
  double  x[6][nwmax], y[6][nwmax], z[6][nwmax];
  // Coordinates for new refinement nodes
  double xn[6][nwmax],yn[6][nwmax],zn[6][nwmax];
  // Gauss-Legendre mapped coordinates
  double xl[21],yl[21],zl[7];

  // ptr[i][j]: ith component (x,y,z) of coordinate for j = n*k for nth node (1..6) for kth element
  ptr = new double* [3];
  for(i = 0; i < 3; i++) ptr[i] = new double[nwmax];
  // ntr[j][i] is the jth node in the ith element
  ntr = new int* [6];
  for(i = 0; i < 6; i++) ntr[i] = new int[nwmax];
  // The dual graph
  //   nne[0][i] is the number of elements touching  global node i
  //   nne[j][i] for j = 1,...,nne[0][i] are the corresponding element labels
  nne = new int* [8];
  for(i = 0; i < 8; i++) nne[i] = new int[nwmax];

  // Create octahedron
  if (nge == 5) {
    x[0][0] = 0;
    y[0][0] = 0;
    z[0][0] = rb;
    x[1][0] = rb;
    y[1][0] = 0;
    z[1][0] = 0;
    x[2][0] = 0;
    y[2][0] = rb;
    z[2][0] = 0;

    x[0][4] = rb;
    y[0][4] = 0;
    z[0][4] = 0;
    x[1][4] = 0;
    y[1][4] = 0;
    z[1][4] = -rb;
    x[2][4] = 0;
    y[2][4] = rb;
    z[2][4] = 0;

    x[0][5] = 0;
    y[0][5] = 0;
    z[0][5] = -rb;
    x[1][5] = -rb;
    y[1][5] = 0;
    z[1][5] = 0;
    x[2][5] = 0;
    y[2][5] = rb;
    z[2][5] = 0;

    x[0][1] = -rb;
    y[0][1] = 0;
    z[0][1] = 0;
    x[1][1] = 0;
    y[1][1] = 0;
    z[1][1] = rb;
    x[2][1] = 0;
    y[2][1] = rb;
    z[2][1] = 0;

    x[0][3] = 0;
    y[0][3] = 0;
    z[0][3] = rb;
    x[1][3] = 0;
    y[1][3] = -rb;
    z[1][3] = 0;
    x[2][3] = rb;
    y[2][3] = 0;
    z[2][3] = 0;

    x[0][7] = rb;
    y[0][7] = 0;
    z[0][7] = 0;
    x[1][7] = 0;
    y[1][7] = -rb;
    z[1][7] = 0;
    x[2][7] = 0;
    y[2][7] = 0;
    z[2][7] = -rb;

    x[0][6] = 0;
    y[0][6] = 0;
    z[0][6] = -rb;
    x[1][6] = 0;
    y[1][6] = -rb;
    z[1][6] = 0;
    x[2][6] = -rb;
    y[2][6] = 0;
    z[2][6] = 0;

    x[0][2] = -rb;
    y[0][2] = 0;
    z[0][2] = 0;
    x[1][2] = 0;
    y[1][2] = -rb;
    z[1][2] = 0;
    x[2][2] = 0;
    y[2][2] = 0;
    z[2][2] = rb;
    nwe = 8;
  } else if( nge == 6 ) {

// Create tandem bridge model

    nwe = 0;
    for(k = 0; k < 2; k++) {
      xl[0] = -rb+k*0.18;
      yl[0] = 0;
      xl[1] = -rb/sqrt(2.0)+k*0.18;
      yl[1] = -rb/sqrt(2.0);
      for( i=2; i<9; i++ ) {
        xl[i] = -0.03+i*0.015+k*0.18;
        yl[i] = -rb;
      }
      xl[9] = 0.09+rb/sqrt(2.0)+k*0.18;
      yl[9] = -rb/sqrt(2.0);
      xl[10] = 0.09+rb+k*0.18;
      yl[10] = 0;
      xl[11] = 0.09+rb/sqrt(2.0)+k*0.18;
      yl[11] = rb/sqrt(2.0);
      for( i=12; i<19; i++ ) {
        xl[i] = 0.27-i*0.015+k*0.18;
        yl[i] = rb;
      }
      xl[19] = -rb/sqrt(2.0)+k*0.18;
      yl[19] = rb/sqrt(2.0);
      xl[20] = -rb+k*0.18;
      yl[20] = 0;
      for( i=0; i<7; i++ ) {
        zl[i] = -0.045+i*0.015;
      }

      for( i=0; i<20; i++ ) {
        for( j=0; j< 6; j++ ) {
          x[0][nwe] = xl[i];
          y[0][nwe] = yl[i];
          z[0][nwe] = zl[j];
          x[1][nwe] = xl[i+1];
          y[1][nwe] = yl[i+1];
          z[1][nwe] = zl[j];
          x[2][nwe] = xl[i];
          y[2][nwe] = yl[i];
          z[2][nwe] = zl[j+1];
          nwe++;
          x[0][nwe] = xl[i+1];
          y[0][nwe] = yl[i+1];
          z[0][nwe] = zl[j+1];
          x[1][nwe] = xl[i];
          y[1][nwe] = yl[i];
          z[1][nwe] = zl[j+1];
          x[2][nwe] = xl[i+1];
          y[2][nwe] = yl[i+1];
          z[2][nwe] = zl[j];
          nwe++;
        }
      }

      n1[0] = 0;
      n1[1] = 1;
      n1[2] = 18;
      n1[3] = 19;
      n2[0] = 8;
      n2[1] = 9;
      n2[2] = 10;
      n2[3] = 11;
      i2[0] = 2;
      i2[1] = 1;
      i3[0] = 1;
      i3[1] = 2;

      for( i=0; i<2; i++ ) {
        iz = 6*i-6;

        for( j=0; j<4; j++ ) {
          x[0][nwe] = xl[n1[j]];
          y[0][nwe] = yl[n1[j]];
          z[0][nwe] = zl[iz];
          x[i2[i]][nwe] = xl[n1[j]+1];
          y[i2[i]][nwe] = yl[n1[j]+1];
          z[i2[i]][nwe] = zl[iz];
          x[i3[i]][nwe] = k*0.18;
          y[i3[i]][nwe] = 0;
          z[i3[i]][nwe] = zl[iz];
          nwe++;
        }

        for( j=0; j<4; j++ ) {
          x[0][nwe] = xl[n2[j]];
          y[0][nwe] = yl[n2[j]];
          z[0][nwe] = zl[iz];
          x[i2[i]][nwe] = xl[n2[j]+1];
          y[i2[i]][nwe] = yl[n2[j]+1];
          z[i2[i]][nwe] = zl[iz];
          x[i3[i]][nwe] = 0.09+k*0.18;
          y[i3[i]][nwe] = 0;
          z[i3[i]][nwe] = zl[iz];
          nwe++;
        }

        for( j=2; j<8; j++ ) {
          x[0][nwe] = xl[j];
          y[0][nwe] = yl[j];
          z[0][nwe] = zl[iz];
          x[i2[i]][nwe] = xl[j+1];
          y[i2[i]][nwe] = yl[j+1];
          z[i2[i]][nwe] = zl[iz];
          x[i3[i]][nwe] = xl[j];
          y[i3[i]][nwe] = 0;
          z[i3[i]][nwe] = zl[iz];
          nwe++;
 
          x[0][nwe] = xl[j+1];
          y[0][nwe] = 0;
          z[0][nwe] = zl[iz];
          x[i2[i]][nwe] = xl[j];
          y[i2[i]][nwe] = 0;
          z[i2[i]][nwe] = zl[iz];
          x[i3[i]][nwe] = xl[j+1];
          y[i3[i]][nwe] = yl[j+1];
          z[i3[i]][nwe] = zl[iz];
          nwe++;
        }

        for( j=12; j<18; j++ ) {
          x[0][nwe] = xl[j];
          y[0][nwe] = yl[j];
          z[0][nwe] = zl[iz];
          x[i2[i]][nwe] = xl[j+1];
          y[i2[i]][nwe] = yl[j+1];
          z[i2[i]][nwe] = zl[iz];
          x[i3[i]][nwe] = xl[j];
          y[i3[i]][nwe] = 0;
          z[i3[i]][nwe] = zl[iz];
          nwe++;
 
          x[0][nwe] = xl[j+1];
          y[0][nwe] = 0;
          z[0][nwe] = zl[iz];
          x[i2[i]][nwe] = xl[j];
          y[i2[i]][nwe] = 0;
          z[i2[i]][nwe] = zl[iz];
          x[i3[i]][nwe] = xl[j+1];
          y[i3[i]][nwe] = yl[j+1];
          z[i3[i]][nwe] = zl[iz];
          nwe++;
        }

      }
    }
  }

  nen = 6;

  // Compute the midpoints
  for( i=0; i<nwe; i++ ) {
    x[3][i] = 0.5*(x[0][i]+x[1][i]);
    y[3][i] = 0.5*(y[0][i]+y[1][i]);
    z[3][i] = 0.5*(z[0][i]+z[1][i]);
    x[4][i] = 0.5*(x[1][i]+x[2][i]);
    y[4][i] = 0.5*(y[1][i]+y[2][i]);
    z[4][i] = 0.5*(z[1][i]+z[2][i]);
    x[5][i] = 0.5*(x[2][i]+x[0][i]);
    y[5][i] = 0.5*(y[2][i]+y[0][i]);
    z[5][i] = 0.5*(z[2][i]+z[0][i]);
  }

  // Fit geometry
  if( nge == 5 ) {
    for( k=0; k<nwe; k++ ) {
      for( l=0; l<6; l++ ) {
        // radius to node
        rad = sqrt(x[l][k]*x[l][k]+y[l][k]*y[l][k]+z[l][k]*z[l][k]);
        // Scale coordinates to [-rb, rb] box
        x[l][k] *= rb/rad;
        y[l][k] *= rb/rad;
        z[l][k] *= rb/rad;
      }
    }
  } else if( nge == 6 ) {
    radc = cos(pi*0.25)*rb-1e-5;
    for( k=0; k<nwe; k++ ) {
      for( l=0; l<6; l++ ) {
        if( x[l][k] < 0 ) {
          rad = sqrt(x[l][k]*x[l][k]+y[l][k]*y[l][k]);
          if( rad > radc ) {
            x[l][k] *= rb/rad;
            y[l][k] *= rb/rad;
          }
        } else if( 0.09 < x[l][k] && x[l][k] < 0.13 ) {
          rad = sqrt((x[l][k]-0.09)*(x[l][k]-0.09)+y[l][k]*y[l][k]);
          if( rad > radc ) {
            x[l][k] = (x[l][k]-0.09)*rb/rad+0.09;
            y[l][k] *= rb/rad;
          }
        } else if( 0.14 < x[l][k] && x[l][k] < 0.18 ) {
          rad = sqrt((x[l][k]-0.18)*(x[l][k]-0.18)+y[l][k]*y[l][k]);
          if( rad > radc ) {
            x[l][k] = (x[l][k]-0.18)*rb/rad+0.18;
            y[l][k] *= rb/rad;
          }
        } else if( 0.27 < x[l][k] ) {
          rad = sqrt((x[l][k]-0.27)*(x[l][k]-0.27)+y[l][k]*y[l][k]);
          if( rad > radc ) {
            x[l][k] = (x[l][k]-0.27)*rb/rad+0.27;
            y[l][k] *= rb/rad;
          }
        }
      }
    }
  }
  if (ndv != 0) {
    // Refine mesh
    for(i = 0; i < ndv; i++) {
      // Number of new elements
      nm = 0;
      for(j = 0; j < nwe; j++) {
        xn[0][nm] = x[0][j];
        yn[0][nm] = y[0][j];
        zn[0][nm] = z[0][j];
        xn[1][nm] = x[3][j];
        yn[1][nm] = y[3][j];
        zn[1][nm] = z[3][j];
        xn[2][nm] = x[5][j];
        yn[2][nm] = y[5][j];
        zn[2][nm] = z[5][j];
        xn[3][nm] = 0.5*(xn[0][nm]+xn[1][nm]);
        yn[3][nm] = 0.5*(yn[0][nm]+yn[1][nm]);
        zn[3][nm] = 0.5*(zn[0][nm]+zn[1][nm]);
        xn[4][nm] = 0.5*(xn[1][nm]+xn[2][nm]);
        yn[4][nm] = 0.5*(yn[1][nm]+yn[2][nm]);
        zn[4][nm] = 0.5*(zn[1][nm]+zn[2][nm]);
        xn[5][nm] = 0.5*(xn[2][nm]+xn[0][nm]);
        yn[5][nm] = 0.5*(yn[2][nm]+yn[0][nm]);
        zn[5][nm] = 0.5*(zn[2][nm]+zn[0][nm]);
        nm++;

        xn[0][nm] = x[3][j];
        yn[0][nm] = y[3][j];
        zn[0][nm] = z[3][j];
        xn[1][nm] = x[1][j];
        yn[1][nm] = y[1][j];
        zn[1][nm] = z[1][j];
        xn[2][nm] = x[4][j];
        yn[2][nm] = y[4][j];
        zn[2][nm] = z[4][j];
        xn[3][nm] = 0.5*(xn[0][nm]+xn[1][nm]);
        yn[3][nm] = 0.5*(yn[0][nm]+yn[1][nm]);
        zn[3][nm] = 0.5*(zn[0][nm]+zn[1][nm]);
        xn[4][nm] = 0.5*(xn[1][nm]+xn[2][nm]);
        yn[4][nm] = 0.5*(yn[1][nm]+yn[2][nm]);
        zn[4][nm] = 0.5*(zn[1][nm]+zn[2][nm]);
        xn[5][nm] = 0.5*(xn[2][nm]+xn[0][nm]);
        yn[5][nm] = 0.5*(yn[2][nm]+yn[0][nm]);
        zn[5][nm] = 0.5*(zn[2][nm]+zn[0][nm]);
        nm++;

        xn[0][nm] = x[5][j];
        yn[0][nm] = y[5][j];
        zn[0][nm] = z[5][j];
        xn[1][nm] = x[4][j];
        yn[1][nm] = y[4][j];
        zn[1][nm] = z[4][j];
        xn[2][nm] = x[2][j];
        yn[2][nm] = y[2][j];
        zn[2][nm] = z[2][j];
        xn[3][nm] = 0.5*(xn[0][nm]+xn[1][nm]);
        yn[3][nm] = 0.5*(yn[0][nm]+yn[1][nm]);
        zn[3][nm] = 0.5*(zn[0][nm]+zn[1][nm]);
        xn[4][nm] = 0.5*(xn[1][nm]+xn[2][nm]);
        yn[4][nm] = 0.5*(yn[1][nm]+yn[2][nm]);
        zn[4][nm] = 0.5*(zn[1][nm]+zn[2][nm]);
        xn[5][nm] = 0.5*(xn[2][nm]+xn[0][nm]);
        yn[5][nm] = 0.5*(yn[2][nm]+yn[0][nm]);
        zn[5][nm] = 0.5*(zn[2][nm]+zn[0][nm]);
        nm++;

        xn[0][nm] = x[3][j];
        yn[0][nm] = y[3][j];
        zn[0][nm] = z[3][j];
        xn[1][nm] = x[4][j];
        yn[1][nm] = y[4][j];
        zn[1][nm] = z[4][j];
        xn[2][nm] = x[5][j];
        yn[2][nm] = y[5][j];
        zn[2][nm] = z[5][j];
        xn[3][nm] = 0.5*(xn[0][nm]+xn[1][nm]);
        yn[3][nm] = 0.5*(yn[0][nm]+yn[1][nm]);
        zn[3][nm] = 0.5*(zn[0][nm]+zn[1][nm]);
        xn[4][nm] = 0.5*(xn[1][nm]+xn[2][nm]);
        yn[4][nm] = 0.5*(yn[1][nm]+yn[2][nm]);
        zn[4][nm] = 0.5*(zn[1][nm]+zn[2][nm]);
        xn[5][nm] = 0.5*(xn[2][nm]+xn[0][nm]);
        yn[5][nm] = 0.5*(yn[2][nm]+yn[0][nm]);
        zn[5][nm] = 0.5*(zn[2][nm]+zn[0][nm]);
        nm++;
      }
      nwe *= 4;
      radc = cos(pi*0.25/i)*rb-1e-5;

      // Copy new points into storage
      for(k = 0; k < nwe; k++) {
        for(l = 0; l < 6; l++) {
          x[l][k] = xn[l][k];
          y[l][k] = yn[l][k];
          z[l][k] = zn[l][k];

          // Fit geometry
          if (nge == 5) {
            rad = sqrt(x[l][k]*x[l][k]+y[l][k]*y[l][k]+z[l][k]*z[l][k]);
            x[l][k] *= rb/rad;
            y[l][k] *= rb/rad;
            z[l][k] *= rb/rad;
          } else if (nge == 6) {
            if( x[l][k] < 0 ) {
              rad = sqrt(x[l][k]*x[l][k]+y[l][k]*y[l][k]);
              if( rad > radc ) {
                x[l][k] *= rb/rad;
                y[l][k] *= rb/rad;
              }
            } else if( 0.09 < x[l][k] && x[l][k] < 0.13 ) {
              rad = sqrt((x[l][k]-0.09)*(x[l][k]-0.09)+y[l][k]*y[l][k]);
              if( rad > radc ) {
                x[l][k] = (x[l][k]-0.09)*rb/rad+0.09;
                y[l][k] *= rb/rad;
              }
            } else if( 0.14 < x[l][k] && x[l][k] < 0.18 ) {
              rad = sqrt((x[l][k]-0.18)*(x[l][k]-0.18)+y[l][k]*y[l][k]);
              if( rad > radc ) {
                x[l][k] = (x[l][k]-0.18)*rb/rad+0.18;
                y[l][k] *= rb/rad;
              }
            } else if( 0.27 < x[l][k] ) {
              rad = sqrt((x[l][k]-0.27)*(x[l][k]-0.27)+y[l][k]*y[l][k]);
              if( rad > radc ) {
                x[l][k] = (x[l][k]-0.27)*rb/rad+0.27;
                y[l][k] *= rb/rad;
              }
            }
          }
          xn[l][k] = 0;
          yn[l][k] = 0;
          zn[l][k] = 0;
        }
      }
    }
  }

  // ptr[j][i]: Coordinate of the jth component for ith global node
  for(i = 0; i < nen; i++) {
    ptr[0][i] = x[i][0];
    ptr[1][i] = y[i][0];
    ptr[2][i] = z[i][0];
    ntr[i][0] = i;
  }
  nwn = nen;

  // ntr[j][i]: Global label on jth node on ith element
  for(i = 1; i < nwe; i++) {
    for(j = 0; j < nen; j++) {
      ic = 0;
      for(k = 0; k < nwn; k++) {
        if (std::abs(x[j][i]-ptr[0][k]) <= eps) {
          if (std::abs(y[j][i]-ptr[1][k]) <= eps) {
            if (std::abs(z[j][i]-ptr[2][k]) <= eps) {
              ic = 1;
              ntr[j][i] = k;
            }
          }
        }
      }
      if (ic == 0) {
        ptr[0][nwn] = x[j][i];
        ptr[1][nwn] = y[j][i];
        ptr[2][nwn] = z[j][i];
        ntr[j][i] = nwn;
        nwn++;
      }
    }
  }
  
  // nne[0][i] is the number of elements touching  global node i
  // nne[j][i] for j = 2,...,nne[0][i]+1 are the corresponding element labels
  for(j = 0; j < 8; j++) {
    for(i = 0; i < nwn; i++) {
      nne[j][i] = 0;
    }
  }
  for(i = 0; i < nwn; i++) {
    nne[0][i] = 0;
    ic = 0;
    for(j = 0; j < nwe; j++) {
      for(k = 0; k < nen; k++) {
        if (ntr[k][j] == i) {
          nne[0][i]++;
          ic++;
          nne[ic][i] = j;
        }
      }
    }
  }
}
