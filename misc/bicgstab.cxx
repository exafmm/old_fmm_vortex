#include "../misc/constants.h"

void matvec(int n, double* x, double (*a)[nwmax], double* b, int it) {
  int i,j;
  double bb;
  for( i = 0; i < n; ++i) {
    bb = 0.0;
    for( j = 0; j < n; ++j) {
      bb += a[i][j]*x[j];
    }
    b[i] = bb;
  }
}

void bicgstab(double* x, double (*a)[nwmax], double* b, int n) {
  bool okprint,GoOn,rcmp,xpdt;
  const char *typestop;
  int i,j,k,info,nmv;
  const double zero = 0.0, one = 1.0, delta = 0.01;
  double kappa0,kappal,maxval1,mxnrmr,mxnrmx;
  double sum1,rnrm0,rnrm,rnrmMax,rnrmin,alpha,omega,sigma,rho0,rho1,beta;
  double varrho,hatgamma;
  double work[10][nwmax],rwork[9][3];

  /*=================
  --- Set options ---
  =================*/

  okprint = false;
  typestop = "rel";

  /*================
  --- Initialize ---
  ================*/

  matvec(n,x,a,work[1],0);
  for( i=0; i<n; i++ ) work[1][i] = b[i]-work[1][i];
  sum1 = zero;
  for( i=0; i<n; i++ ) {
    work[0][i] = work[1][i];
    work[8][i] = work[1][i];
    work[7][i] = x[i];
    sum1 += work[1][i]*work[1][i];
  }
  rnrm0 = sqrt(sum1);
  rnrm = rnrm0;
  if (strcmp(typestop, "max")) {
    maxval1 = zero;
    for( i=0; i<n; i++ ) maxval1 = std::max(maxval1,std::abs(work[1][i]));
    rnrmMax = maxval1;
  }
  info = 0;
  nmv = 1;
  rnrmin = rnrm;
  mxnrmx = rnrm0;
  mxnrmr = rnrm0;
  rcmp = false;
  xpdt = false;
  alpha = zero;
  omega = one;
  sigma = one;
  rho0 =  one;

  /*===============================
  --- Check initial convergence ---
  ===============================*/

  if (strcmp(typestop, "rel")) {
    GoOn = rnrm >= tol*rnrm0;
  } else if (strcmp(typestop, "abs")) {
    GoOn = rnrm >= tol;
  } else if (strcmp(typestop, "max")) {
    GoOn = rnrmMax >= tol;
  }

  /*#####################
  --- Begin Main loop ---
  #####################*/

  while( GoOn ) {

  /*===================
  --- The BiCG part ---
  ===================*/

    rho0 = -omega*rho0;
    for( i=0; i<2; i++ ) {
      sum1 = zero;
      for( j=0; j<n; j++ ) sum1 += work[0][j]*work[i+1][j];
      rho1 = sum1;
      if( rho0 == zero ) {
        info = 2;
        return;
      }
      beta = alpha*(rho1/rho0);
      rho0 = rho1;
      for( j=0; j<=i; j++ ) {
        for( k=0; k<n; k++ ) {
          work[j+4][k] = work[j+1][k]-beta*work[j+4][k];
        }
      }
      matvec(n,work[i+4],a,work[i+5],1);
      nmv++;
      sum1 = zero;
      for( j=0; j<n; j++ ) sum1 += work[0][j]*work[i+5][j];
      sigma = sum1;
      if (sigma == 0.0) {
        info = 2;
        return;
      }
      alpha = rho1/sigma;
      for( j=0; j<n; j++ ) x[j] = alpha*work[4][j]+x[j];
      for( j=0; j<=i; j++ ) {
        for( k=0; k<n; k++ ) {
          work[j+1][k] = -alpha*work[j+5][k]+work[j+1][k];
        }
      }
      matvec(n,work[i+1],a,work[i+2],1);
      nmv++;
      sum1 = zero;
      for( j=0; j<n; j++ ) sum1 += work[1][j]*work[1][j];
      rnrm = sqrt(sum1);
      if (strcmp(typestop, "max")) {
        maxval1 = zero;
        for( j=0; j<n; j++ ) maxval1 = std::max(maxval1,std::abs(work[1][j]));
        rnrmMax = maxval1;
      }
      mxnrmx = std::max(mxnrmx,rnrm);
      mxnrmr = std::max(mxnrmr,rnrm);
    }

    /*================================
    --- The convex polynomial part ---
    ================================*/

    for( i=0; i<3; i++ ) {
      for( j=i-1; j<2; j++ ) {
        sum1 = zero;
        for( k=0; k<n; k++ ) {
          sum1 += work[j+2][k]*work[i+1][k];
        }
        rwork[i][j+1] = sum1;
        rwork[j+1][i] = rwork[i][j+1];
      }
    }
    for( i=3; i<6; i++ ) for( j=0; j<3; j++ ) rwork[i][j] = rwork[j][i-3];
    rwork[6][0] = -one;
    rwork[6][1] = rwork[0][1]/rwork[4][1];
    rwork[6][2] = zero;
    rwork[7][0] = zero;
    rwork[7][1] = rwork[2][1]/rwork[4][1];
    rwork[7][2] = -one;
    for( i=0; i<3; i++ ) rwork[8][i] = zero;
    for( i=0; i<3; i++ ) for( j=0; j<3; j++ ) rwork[8][j] += rwork[7][i]*rwork[i][j];
    sum1 = zero;
    for( i=0; i<3; i++ ) sum1 += rwork[7][i]*rwork[8][i];
    kappal = sqrt(sum1);
    for( i=0; i<3; i++ ) rwork[8][i] = zero;
    for( i=0; i<3; i++ ) for( j=0; j<3; j++ ) rwork[8][j] += rwork[6][i]*rwork[i][j];
    sum1 = zero;
    for( i=0; i<3; i++ ) sum1 += rwork[6][i]*rwork[8][i];
    kappa0 = sqrt(sum1);
    sum1 = zero;
    for( i=0; i<3; i++ ) sum1 += rwork[7][i]*rwork[8][i];
    varrho = sum1;
    varrho /= kappa0*kappal;
    hatgamma = (varrho > 0 ? 1 : (varrho < 0 ? -1 : 0))*
               std::max(std::abs(varrho),0.7)*(kappa0/kappal);
    for( i=0; i<3; i++ ) rwork[6][i] -= hatgamma*rwork[7][i];

    /*======================
    --- Update variables ---
    ======================*/

    omega = rwork[6][2];
    for( i=0; i<2; i++ ) {
      for( j=0; j<n; j++ ) {
        work[4][j] -= rwork[6][i+1]*work[i+5][j];
        x[j] += rwork[6][i+1]*work[i+1][j];
        work[1][j] -= rwork[6][i+1]*work[i+2][j];
      }
    }
    for( i=0; i<3; i++ ) rwork[8][i] = zero;
    for( i=0; i<3; i++ ) {
      for ( j=0; j<3; j++ ) {
        rwork[8][j] += rwork[6][i]*rwork[i][j];
      }
    }
    sum1 = zero;
    for( i=0; i<3; i++ ) sum1 += rwork[6][i]*rwork[8][i];
    rnrm = sqrt(sum1);

    /*==============================
    --- The reliable update part ---
    ==============================*/

    mxnrmx = std::max(mxnrmx,rnrm);
    mxnrmr = std::max(mxnrmr,rnrm);
    xpdt = rnrm < delta*rnrm0 && rnrm0 < mxnrmx;
    rcmp = (rnrm < delta*mxnrmr && rnrm0 < mxnrmr) || xpdt;
    if( rcmp ) {
      matvec(n,x,a,work[1],1);
      nmv++;
      for( i=0; i<n; i++ ) work[1][i] = work[8][i]-work[1][i];
      mxnrmr = rnrm;
      if( xpdt ) {
        for( i=0; i<n; i++ ) {
          work[7][i]+= x[i];
          x[i] = zero;
          work[8][i] = work[1][i];
        }
        mxnrmx = rnrm;
      }
    }
    if( rnrm < rnrmin && std::abs(x[0]) < 1 ) {
      rnrmin = rnrm;
      for( i=0; i<n; i++ ) work[9][i] = work[7][i]+x[i];
      if (strcmp(typestop, "rel")) {
        if( okprint ) std::cout << nmv << " " << rnrmin/rnrm0 << std::endl;
      } else if (strcmp(typestop, "abs")) {
        if( okprint ) std::cout << nmv << " " << rnrmin << std::endl;
      } else if (strcmp(typestop, "max")) {
        maxval1 = zero;
        for( i=0; i<n; i++ ) maxval1 = std::max(maxval1,std::abs(work[1][i]));
        rnrmMax = maxval1;
        if( okprint ) std::cout << nmv << " " << rnrmMax << std::endl;
      }
    }

    /*=======================
    --- Check convergence ---
    =======================*/

    if (strcmp(typestop, "rel")) {
      GoOn = rnrmin >= tol*rnrm0 && nmv < itmax;
    } else if (strcmp(typestop, "abs")) {
      GoOn = rnrmin >= tol && nmv < itmax;
    } else if (strcmp(typestop, "max")) {
      maxval1 = zero;
      for( i=0; i<n; i++ ) maxval1 = std::max(maxval1,std::abs(work[1][i]));
      rnrmMax = maxval1;
      GoOn = rnrmMax >= tol && nmv < itmax;
    }

    /*###################
    --- End main loop ---
    ###################*/

    for( i=0; i<n; i++ ) x[i] = work[9][i];
  }

}
