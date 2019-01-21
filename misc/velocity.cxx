#include "../misc/constants.h"

extern float *xi,*yi,*zi,*uw,*vw,*ww;
extern int    **ntr,*nne;
extern double **ptr,**vna,**vta,**vsa,**p1,**p2,**p3,**p4,**p5,**p6,**hsg,*are,*crv,*si,*gt,*gs;

void velocity(int n0, int n1,int nsv, int np, int nge) {
  int l,k,j,i,ii;
  double x0,y0,z0,u0,v0,w0,x,y,z,siint,gxint,gyint,gzint,dxij,dyij,dzij,r,den,Gx,Gy,Gz,cf;
  double ph[6];

  for( l=0; l<np; l++ ) {
    x0 = xi[l];
    y0 = yi[l];
    z0 = zi[l];
    u0 = 0.0;
    v0 = 0.0;
    w0 = 0.0;
    for( k=n0; k<n1; k++ ) {
      for( j=0; j<ngl; j++ ) {
        x = 0.0; 
        y = 0.0; 
        z = 0.0; 
        siint = 0.0;
        gxint = 0.0;
        gyint = 0.0;
        gzint = 0.0;

        ph[0] = p1[j][k];
        ph[1] = p2[j][k];
        ph[2] = p3[j][k];
        ph[3] = p4[j][k];
        ph[4] = p5[j][k];
        ph[5] = p6[j][k];

        for( i=0; i<6; i++ ) {
          ii = ntr[i][k];

          x += ptr[0][ii]*ph[i];
          y += ptr[1][ii]*ph[i];
          z += ptr[2][ii]*ph[i];
          siint += si[ii]*ph[i];
          gxint += gs[ii]*vsa[0][ii]*ph[i]-gt[ii]*vta[0][ii]*ph[i];
          gyint += gs[ii]*vsa[1][ii]*ph[i]-gt[ii]*vta[1][ii]*ph[i];
          gzint += gs[ii]*vsa[2][ii]*ph[i]-gt[ii]*vta[2][ii]*ph[i];
        }

// Compute the Green's function derivative

        dxij = x0-x;
        dyij = y0-y;
        dzij = z0-z;

        r = sqrt(dxij*dxij+dyij*dyij+dzij*dzij);
        den = 4*pi*r*r*r;

        Gx = dxij/den;
        Gy = dyij/den;
        Gz = dzij/den;

// Apply the triangle quadrature

        cf = 0.5*hsg[j][k];

        if( nsv == 1 ) {
          u0 += siint*Gx*cf;
          v0 += siint*Gy*cf;
          w0 += siint*Gz*cf;
        } else if( nsv == 2 ) {
          u0 += (gyint*Gz-gzint*Gy)*cf;
          v0 += (gzint*Gx-gxint*Gz)*cf;
          w0 += (gxint*Gy-gyint*Gx)*cf;
        }
      }
    }
    uw[l] = u0;
    vw[l] = v0;
    ww[l] = w0;
  }
}
