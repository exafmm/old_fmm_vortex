#include "../misc/constants.h"

extern int *nej,*nfj,*nlbj,*nel,*nlbl,*nem,*nlbm,**npx,*njb;

extern void jsbox(int, int, int, int&, int, int, int);
extern void ijbox(int, int, int, int, int);

void isbox(int n0, int n1, int mj, int& mjp, int mp, int lbi, int lbj, int lbl, int lbm, int& lbjr, int& lbjrp, int lev, int ipb, int npb, int jlm, double* tfmm) {
  int i,j,jj,lbjp;
  double tic,toc;

  tic = get_time();

  if( lev!=1 && ipb != -2 ) {
    for( i=0; i<nbnp; i++ ) {
      for( j=0; j<3; j++ ) {
        npx[j][i] = 0;
      }
    }
    jsbox(lbj,lbl,lbm,lbjp,lev,ipb,npb);
    ijbox(lbi,lbjp,lev,ipb,npb);
    lbjr = std::max(lbj,lbl);
    lbjr = std::max(lbjr,lbm);
    lbjrp = lbjp;
    for( jj=0; jj<lbjrp; jj++ ) {
      if( npx[2][jj] == -1 ) {
        njb[jj] = nem[nfj[jj]]+nlbm[lev-1];
      } else if ( npx[2][jj] == 1 ) {
        njb[jj] = nel[nfj[jj]]+nlbl[lev-1];
      } else {
        njb[jj] = nej[nfj[jj]]+nlbj[lev-1];
      }
    }
  }

  toc = tic;
  tic = get_time();
  tfmm[8] += tic-toc;

}
