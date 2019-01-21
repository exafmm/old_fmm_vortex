#include "../misc/constants.h"

extern double **dima;

extern void ludcmp(double**, int, int*);
extern void lubksb(double**, int, int*, double*, int);

void inv(double** a, int n){
  int i,j,indx[nwmax];

  dima = new double* [nwmax];
  for( i=0; i<nwmax; i++ ) dima[i] = new double [nwmax];

  for( i=0; i<n; i++ ) {
    for( j=0; j<n; j++ ) {
      dima[i][j] = 0;
    }
    dima[i][i] = 1;
  }

  ludcmp(a,n,indx);

  for( i=0; i<n; i++ ) {
    lubksb(a,n,indx,dima[i],i);
  }

  for( i=0; i<n; i++ ) {
    for( j=0; j<n; j++ ) {
      a[i][j] = dima[i][j];
    }
  }

  for( i=0; i<nwmax; i++ ) delete[] dima[i];
  delete[] dima;
}

// LU decomposition

void ludcmp(double** a, int n, int* indx) {
  int i,j,k,imax;
  double aamax,asum,dum,vv[nwmax];

  for( i=0; i<n; i++ ) {
    aamax = 0;
    for( j=0; j<n; j++ ) {
      aamax = std::max(std::abs(a[j][i]),aamax);
    }
    if( aamax < eps ) std::cout << "singular matrix in ludcmp" << std::endl;
    vv[i] = 1.0/aamax;
  }

  for( j=0; j<n; j++ ) {

    for( i=0; i<j; i++ ) {
      asum = a[j][i];
      for( k=0; k<i; k++ ) {
        asum -= a[k][i]*a[j][k];
      }
      a[j][i] = asum;
    }

    aamax = 0;
    for( i=j; i<n; i++ ) {
      asum = a[j][i];
      for( k=0; k<j; k++ ) {
        asum -= a[k][i]*a[j][k];
      }
      a[j][i] = asum;
      dum = vv[i]*std::abs(asum);
      if( dum >= aamax ) {
        imax = i;
        aamax = dum;
      }
    }

    if( j != imax ) {
      for( k=0; k<n; k++ ) {
        dum = a[k][imax];
        a[k][imax] = a[k][j];
        a[k][j] = dum;
      }
      vv[imax] = vv[j];
    }

    indx[j] = imax;
    if( a[j][j] == 0 ) a[j][j] = eps;
    if( j != n ) {
      dum = 1.0/a[j][j];
      for( i=j+1; i<n; i++ ) a[j][i] *= dum;
    }
  }
}

void lubksb(double** a, int n, int* indx, double* b, int k) {
  int i,j,ii,ll;
  double asum;

  ii = -1;
  for( i=0; i<n; i++ ) {
    ll = indx[i];
    asum = b[ll];
    b[ll] = b[i];
    if( ii != -1 ) {
      for( j=ii; j<i; j++ ) asum -= a[j][i]*b[j];
    } else {
      ii = i;
    }
    b[i] = asum;
  }

  for( i=n-1; i>=0; i-- ) {
    asum = b[i];
    for( j=i+1; j<n; j++ ) asum -= a[j][i]*b[j];
    b[i] = asum/a[i][i];
  }
}
