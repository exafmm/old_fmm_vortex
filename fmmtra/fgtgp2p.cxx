#include "../misc/constants.h"

extern float *xi,*yi,*zi,*xj,*yj,*zj,*vj,*sj;
extern int **ndi,**ndj,*nej,*nij,**neij;
extern float *gxd,*gyd;

void fgtp2p(int lbi, int lbj) {
  int ii,ij,jj,i,j;
  double gij,dxij,dyij,dzij,rij,sij,cutoff;

  for( ii=0; ii<lbi; ii++ ) {
    for( ij=0; ij<nij[ii]; ij++ ) {
      jj = neij[ij][ii];
      for( i=ndi[0][ii]; i<=ndi[1][ii]; i++ ) {
        gij = 0;
        for( j=ndj[0][jj]; j<=ndj[1][jj]; j++ ) {
          #include "../direct/rij.cxx"
          sij = 2*sj[j]*sj[j];
          if( rij < 100*sij ) {
            cutoff = 1/pow(pi*sij,1.5)*exp(-rij/sij);
            gij += gyd[j]*cutoff;
          }
        }
        gxd[i] += gij;
      } 
    }
  }
}
