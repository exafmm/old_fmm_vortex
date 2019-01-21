#include "../misc/constants.h"

extern int **npx;

extern void ijbox(int, int, int, int, int);

void ifbox(int n0, int n1, int& mjp, int mp, int lbi, int lbj, int& lbjr, int& lbjrp, int lev, int ipb, int npb, double* tfmm) {
  int i,j;
  double tic,toc;

  tic = get_time();

  if( ipb != -2 ) {
    for( i=0; i<nbnp; i++ ) {
      for( j=0; j<3; j++ ) {
        npx[j][i] = 0;
      }
    }
    ijbox(lbi,lbj,lev,ipb,npb);
    lbjr = lbj;
    lbjrp = lbj;
  }

  toc = tic;
  tic = get_time();
  tfmm[8] += tic-toc;

}
