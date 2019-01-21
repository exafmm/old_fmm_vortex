#include "../misc/constants.h"

extern float *xi,*yi,*zi,*gxi,*gyi,*gzi;
extern float *xj,*yj,*zj,*gxj,*gyj,*gzj,*vj,*sj;
extern float *gxd,*gyd,*gzd;

extern void memoryuse();
extern void memoryfree();

void dird(int n0, int n1, int n2, int n3) {
  int i,j;
  double gxdd,gydd,gzdd,dxij,dyij,dzij,rij,sij,rsij,cutoff,st,tr;

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
      sij = 2*sj[j]*sj[j];
      rsij = rij/sij;
      cutoff = (3*erf(sqrt(rsij))-(2*rsij+3)*sqrt(4.0/pi*rsij)*exp(-rsij))/3.0/pow(rij,2.5);
      st = gxi[i]*dxij+gyi[i]*dyij+gzi[i]*dzij;
      tr = gxi[i]*(gyj[j]*dzij-gzj[j]*dyij)+gyi[i]*(gzj[j]*dxij-gxj[j]*dzij)+gzi[i]*(gxj[j]*dyij-gyj[j]*dxij);
      gxdd += 0.375/pi*((gyj[j]*dzij-gzj[j]*dyij)*st+dxij*tr)*cutoff;
      gydd += 0.375/pi*((gzj[j]*dxij-gxj[j]*dzij)*st+dyij*tr)*cutoff;
      gzdd += 0.375/pi*((gxj[j]*dyij-gyj[j]*dxij)*st+dzij*tr)*cutoff;
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
