#include "../misc/constants.h"

extern float *xi,*yi,*zi,*gxi,*gyi,*gzi;
extern float *xj,*yj,*zj,*gxj,*gyj,*gzj,*vj,*sj;
extern float *gxd,*gyd,*gzd;

extern void memoryuse();
extern void memoryfree();

void dir(int n0, int n1, int n2, int n3) {
  int i,j;
  double gxdd,gydd,gzdd,dxij,dyij,dzij,rij,sij,cutoff;

  gxd = new float [npmax];
  gyd = new float [npmax];
  gzd = new float [npmax];
  mem = npmax*3*4;
  memoryuse();

  for( i=n0; i<=n1; i++ ) {
    gxdd = 0;
    gydd = 0;
    gzdd = 0;
    for( j=n2; j<=n3; j++ ){
      #include "rij.cxx"
      sij = sj[j]*sj[j];
      cutoff = (rij+2.5*sij)/pow(rij+sij,2.5);
      gxdd += 0.25/pi*(gyi[i]*gzj[j]-gzi[i]*gyj[j])*cutoff;
      gydd += 0.25/pi*(gzi[i]*gxj[j]-gxi[i]*gzj[j])*cutoff;
      gzdd += 0.25/pi*(gxi[i]*gyj[j]-gyi[i]*gxj[j])*cutoff;
      cutoff = (rij+3.5*sij)/pow(rij+sij,3.5)*(gxi[i]*dxij+gyi[i]*dyij+gzi[i]*dzij);
      gxdd += 0.75/pi*(gyj[j]*dzij-gzj[j]*dyij)*cutoff;
      gydd += 0.75/pi*(gzj[j]*dxij-gxj[j]*dzij)*cutoff;
      gzdd += 0.75/pi*(gxj[j]*dyij-gyj[j]*dxij)*cutoff;
    }
    gxd[i] = gxdd;
    gyd[i] = gydd;
    gzd[i] = gzdd;
  }
  for( i=0; i<npmax; i++ ) {
    gxi[i] = 0;
    gyi[i] = 0;
    gzi[i] = 0;
  }
  for( i=n0; i<=n1; i++ ) {
    gxi[i] = gxd[i];
    gyi[i] = gyd[i];
    gzi[i] = gzd[i];
  }

  delete[] gxd;
  delete[] gyd;
  delete[] gzd;
  mem = npmax*3*4;
  memoryfree();
}
