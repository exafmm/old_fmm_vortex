#include "../misc/constants.h"

extern float *xi,*yi,*zi,*gxi,*gyi,*gzi;
extern float *xj,*yj,*zj,*gxj,*gyj,*gzj,*vj,*sj;
extern float *gxd,*gyd;

extern void memoryuse();
extern void memoryfree();

void dird(int n0, int n1, int n2, int n3) {
  int i,j;
  double gd,dxij,dyij,dzij,rij,sij,cutoff;

  gxd = new float [npmax];
  gyd = new float [npmax];
  mem = npmax*2*4;
  memoryuse();

  if( neq==5 ) {
    for( i=n2; i<=n3; i++ ) gyd[i] = gxj[i];
  } else if( neq==6 ) {
    for( i=n2; i<=n3; i++ ) gyd[i] = gyj[i];
  } else if( neq==7 ) {
    for( i=n2; i<=n3; i++ ) gyd[i] = gzj[i];
  }

  for( i=n0; i<=n1; i++ ) {
    gd = 0;
    for( j=n2; j<=n3; j++ ) {
      #include "../direct/rij.cxx"
      sij = 2*sj[j]*sj[j];
      cutoff = 1/pow(pi*sij,1.5)*exp(-rij/sij);
      gd += gyd[j]*cutoff;
    }
    gxd[i] = gd;
  }

  if( neq==5 ) {
    for( i=0; i<npmax; i++ ) gxi[i] = 0;
    for( i=n0; i<=n1; i++ ) gxi[i] = gxd[i];
  } else if( neq==6 ) {
    for( i=0; i<npmax; i++ ) gyi[i] = 0;
    for( i=n0; i<=n1; i++ ) gyi[i] = gxd[i];
  } else if( neq==7 ) {
    for( i=0; i<npmax; i++ ) gzi[i] = 0;
    for( i=n0; i<=n1; i++ ) gzi[i] = gxd[i];
  }

  delete[] gxd;
  delete[] gyd;
  mem = npmax*2*4;
  memoryfree();
}
