#include "../misc/parameters.h"
#include "../misc/constants.h"

extern float *xi,*yi,*zi,*gxi,*gyi,*gzi,*vi;
extern float *xj,*yj,*zj,*gxj,*gyj,*gzj,*vj,*sj;
extern float *gxd,*gyd,*gzd;

extern void memoryuse();
extern void memoryfree();

void dir(int n0, int n1, int n2, int n3) {
  int i,j;
  double gxdd,gydd,gzdd,dxij,dyij,dzij,rij,sij,rsij,cutoff;

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
      cutoff = 2.0*vis*pow(sij,2.0)*105/8.0/pi/pow(rij+sij,4.5);
      gxdd +=(gxj[j]*vi[i]-gxi[i]*vj[j])*cutoff;
      gydd +=(gyj[j]*vi[i]-gyi[i]*vj[j])*cutoff;
      gzdd +=(gzj[j]*vi[i]-gzi[i]*vj[j])*cutoff;
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
