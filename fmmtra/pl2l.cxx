#include "../misc/constants.h"

extern double **psm;

void pl2l(double* ppx, double* ppy, double* ppz, double* ppxd, double* ppyd, double* ppzd, int nmp, int ni) {
  int j,k,njk;
  double ppxdd,ppydd,ppzdd;

  for( j=0; j<nmp; j++ ) {
    ppxdd = 0;
    ppydd = 0;
    ppzdd = 0;
    for( k=0; k<nmp; k++ ) {
      njk = j*nmp+k;
      ppxdd += ppxd[k]*psm[ni][njk];
      ppydd += ppyd[k]*psm[ni][njk];
      ppzdd += ppzd[k]*psm[ni][njk];
    }
    ppx[j] = ppxdd;
    ppy[j] = ppydd;
    ppz[j] = ppzdd;
  }
}
