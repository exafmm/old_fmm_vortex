#include "../misc/constants.h"

extern double **ptr;
extern int **ntr;
extern double **vna,**vta,**vsa,**p1,**p2,**p3,**p4,**p5,**p6,**hsg;

void matvec2(int nwe, int nwn, int ngl, int nsv,  double* xx, double** dimh, double* bb) {
  int i,m,mc,l,k,j,i1,i2,i3,i4,i5,i6,ii;
  double x0,y0,z0,q0,ptl,ptl1,ptl2,ptl3,ptl4,test;
  double x,y,z,vnx,vny,vnz,vtx,vty,vtz,vsx,vsy,vsz,qint;
  double dxij,dyij,dzij,r,den,cf,Gn,Gtt,Gts,Gst,Gss,bbb;
  double q[nwmax],ph[6];

  for( i=0; i<nwn; i++ ) {
    dimh[i][i] -= 0.5;
  }
  if( nsv == 2 ) {
    for( i=nwn; i<2*nwn; i++ ) {
      dimh[i][i] -= 0.5;
    }
  }
  for( i=0; i<nwn; i++ ) q[i] = 0;
  for( m=0; m<nwn; m++ ) {          // Run over source nodes
    mc = (int) ((double) m)/nwn*100;
    q[m] = 1;                       // Impulse
    for( l=0; l<nwn; l++ ) {        // Run over target nodes
      x0 = ptr[0][l];
      y0 = ptr[1][l];
      z0 = ptr[2][l];
      q0 = q[l];
      ptl = 0;
      ptl1 = 0;
      ptl2 = 0;
      ptl3 = 0;
      ptl4 = 0;
      for( k=0; k<nwe; k++ ) {      // Run over source elements
        i1 = ntr[0][k];
        i2 = ntr[1][k];
        i3 = ntr[2][k];
        i4 = ntr[3][k];
        i5 = ntr[4][k];
        i6 = ntr[5][k];
        test = std::abs(q[i1])+std::abs(q[i2])+std::abs(q[i3])+std::abs(q[i4])
              +std::abs(q[i5])+std::abs(q[i6])+std::abs(q0);
        if( test > tol ) {
          for( j=0; j<ngl; j++ ) {  // Run over quadratures
            x    = 0;
            y    = 0;
            z    = 0;
            vnx  = 0;
            vny  = 0;
            vnz  = 0;
            vtx  = 0;
            vty  = 0;
            vtz  = 0;
            vsx  = 0;
            vsy  = 0;
            vsz  = 0;
            qint = 0;

            ph[0] = p1[j][k];
            ph[1] = p2[j][k];
            ph[2] = p3[j][k];
            ph[3] = p4[j][k];
            ph[4] = p5[j][k];
            ph[5] = p6[j][k];
            for( i=0; i<6; i++ ) {  // Run over node contributions
              ii = ntr[i][k];
              x += ptr[0][ii]*ph[i];
              y += ptr[1][ii]*ph[i];
              z += ptr[2][ii]*ph[i];
              vnx += vna[0][ii]*ph[i];
              vny += vna[1][ii]*ph[i];
              vnz += vna[2][ii]*ph[i];
              vtx += vta[0][ii]*ph[i];
              vty += vta[1][ii]*ph[i];
              vtz += vta[2][ii]*ph[i];
              vsx += vsa[0][ii]*ph[i];
              vsy += vsa[1][ii]*ph[i];
              vsz += vsa[2][ii]*ph[i];
              qint += q[ii]*ph[i];
            }

// Compute the Green's function derivative and apply the triangle quadrature

            dxij = x0-x;
            dyij = y0-y;
            dzij = z0-z;
            r = sqrt(dxij*dxij+dyij*dyij+dzij*dzij);
            den = 4*pi*r*r*r;
            cf = 0.5*hsg[j][k];

            if( nsv == 1 ) {
              Gn = vnx*dxij+vny*dyij+vnz*dzij;
              ptl += (qint-q0)*Gn/den*cf;
            } else if( nsv == 2 ) {
              Gtt = vta[0][l]*(vty*dzij-vtz*dyij)+vta[1][l]*(vtz*dxij-vtx*dzij)+vta[2][l]*(vtx*dyij-vty*dxij);
              Gts = vta[0][l]*(vsy*dzij-vsz*dyij)+vta[1][l]*(vsz*dxij-vsx*dzij)+vta[2][l]*(vsx*dyij-vsy*dxij);
              Gst = vsa[0][l]*(vty*dzij-vtz*dyij)+vsa[1][l]*(vtz*dxij-vtx*dzij)+vsa[2][l]*(vtx*dyij-vty*dxij);
              Gss = vsa[0][l]*(vsy*dzij-vsz*dyij)+vsa[1][l]*(vsz*dxij-vsx*dzij)+vsa[2][l]*(vsx*dyij-vsy*dxij);
              ptl1 += (qint-q0)*Gtt/den*cf;
              ptl2 += (qint-q0)*Gts/den*cf;
              ptl3 += (qint-q0)*Gst/den*cf;
              ptl4 += (qint-q0)*Gss/den*cf;
            }

          }
        }
      }
      if( nsv == 1 ) {
        dimh[m][l] = ptl-0.5*q0;
      } else if( nsv == 2 ) {
        dimh[m][l] = ptl3-0.5*q0;
        dimh[m+nwn][l] = -ptl4+0.5*q0;
        dimh[m][l+nwn] = ptl1-0.5*q0;
        dimh[m+nwn][l+nwn] = -ptl2+0.5*q0;
      }
    }
    q[m] = 0; // Reset
  }

  for( i=0; i<nwn; i++ ) {
    bbb = 0;
    for( j=0; j<nwn; j++ ) {
      bbb += dimh[i][j]*xx[j];
    }
    bb[i] = bbb;
  }
}

void bicgstab2(double x[], double** a, double b[], int nwe, int nwn, int ngl, int nsv) {
  bool okprint,GoOn,rcmp,xpdt;
  const char *typestop;
  int i,j,k,info,nmv;
  double zero = 0, one = 1, delta = 0.01;
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

  matvec2(nwe,nwn,ngl,nsv,x,a,work[1]);
  for( i=0; i<nwn; i++ ) work[1][i] = b[i]-work[1][i];
  sum1 = zero;
  for( i=0; i<nwn; i++ ) {
    work[0][i] = work[1][i];
    work[8][i] = work[1][i];
    work[7][i] = x[i];
    sum1 += work[1][i]*work[1][i];
  }
  rnrm0 = sqrt(sum1);
  rnrm = rnrm0;
  if (strcmp(typestop, "max")) {
    maxval1 = zero;
    for( i=0; i<nwn; i++ ) maxval1 = std::max(maxval1,std::abs(work[1][i]));
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
      for( j=0; j<nwn; j++ ) sum1 += work[0][j]*work[i+1][j];
      rho1 = sum1;
      if( rho0 == zero ) {
        info = 2;
        return;
      }
      beta = alpha*(rho1/rho0);
      rho0 = rho1;
      for( j=0; j<=i; j++ ) {
        for( k=0; k<nwn; k++ ) {
          work[j+4][k] = work[j+1][k]-beta*work[j+4][k];
        }
      }
      matvec2(nwe,nwn,ngl,nsv,work[i+4],a,work[i+5]);
      nmv++;
      sum1 = zero;
      for( j=0; j<nwn; j++ ) sum1 += work[0][j]*work[i+5][j];
      sigma = sum1;
      if( sigma == zero ) {
        info = 2;
        return;
      }
      alpha = rho1/sigma;
      for( j=0; j<nwn; j++ ) x[j] = alpha*work[4][j]+x[j];
      for( j=0; j<=i; j++ ) {
        for( k=0; k<nwn; k++ ) {
          work[j+1][k] = -alpha*work[j+5][k]+work[j+1][k];
        }
      }
      matvec2(nwe,nwn,ngl,nsv,work[i+1],a,work[i+2]);
      nmv++;
      sum1 = zero;
      for( j=0; j<nwn; j++ ) sum1 += work[1][j]*work[1][j];
      rnrm = sqrt(sum1);
      if (strcmp(typestop, "max")) {
        maxval1 = zero;
        for( i=0; i<nwn; i++ ) maxval1 = std::max(maxval1,std::abs(work[1][i]));
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
        for( k=0; k<nwn; k++ ) {
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
      for( j=0; j<nwn; j++ ) {
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
      matvec2(nwe,nwn,ngl,nsv,x,a,work[1]);
      nmv++;
      for( i=0; i<nwn; i++ ) work[1][i] = work[8][i]-work[1][i];
      mxnrmr = rnrm;
      if( xpdt ) {
        for( i=0; i<nwn; i++ ) {
          work[7][i]+= x[i];
          x[i] = zero;
          work[8][i] = work[1][i];
        }
        mxnrmx = rnrm;
      }
    }
    if( rnrm < rnrmin && std::abs(x[0]) < 1 ) {
      rnrmin = rnrm;
      for( i=0; i<nwn; i++ ) work[9][i] = work[7][i]+x[i];
      if (strcmp(typestop, "rel")) {
        if( okprint ) std::cout << nmv << " " << rnrmin/rnrm0 << std::endl;
      } else if (strcmp(typestop, "abs")) {
        if( okprint ) std::cout << nmv << " " << rnrmin << std::endl;
      } else if (strcmp(typestop, "max")) {
        maxval1 = zero;
        for( i=0; i<nwn; i++ ) maxval1 = std::max(maxval1,std::abs(work[1][i]));
        rnrmMax = maxval1;
        if( okprint ) std::cout << nmv << " " << rnrmMax << std::endl;
      }
    }

    /*=======================
    --- Check convergence ---
    =======================*/

    if (strcmp(typestop, "rel")) {
      GoOn = rnrm >= tol*rnrm0 && nmv < itmax;
    } else if (strcmp(typestop, "abs")) {
      GoOn = rnrm >= tol && nmv < itmax;
    } else if (strcmp(typestop, "max")) {
      maxval1 = zero;
      for( i=0; i<nwn; i++ ) maxval1 = std::max(maxval1,std::abs(work[1][i]));
      rnrmMax = maxval1;
      GoOn = rnrmMax >= tol && nmv < itmax;
    }

    /*###################
    --- End main loop ---
    ###################*/

    for( i=0; i<nwn; i++ ) x[i] = work[9][i];
  }

}
