#include "../misc/constants.h"

extern int **npx;

extern void jpbox(int, int&, int, int, int);
extern void ijbox(int, int, int, int, int);

void ipbox(int n0, int n1, int& mjp, int mp, int lbi, int lbj, int& lbjr, int& lbjrp, int lev, int ipb, int npb, double* tfmm) {
  int i,j,lbjp;
  double tic,toc;

  tic = get_time();

  if( lev!=1 && ipb != -2 ) {
    for( i=0; i<nbnp; i++ ) {
      for( j=0; j<3; j++ ) {
        npx[j][i] = 0;
      }
    }
    jpbox(lbj,lbjp,lev,ipb,npb);
    ijbox(lbi,lbjp,lev,ipb,npb);
    lbjr = lbj;
    lbjrp = lbjp;
  }

  toc = tic;
  tic = get_time();
  tfmm[8] += tic-toc;

}
