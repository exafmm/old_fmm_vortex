#include "mpi.h"
#include "../misc/constants.h"

extern float *xj,*yj,*zj,*gxj,*gyj,*gzj,*vj,*sj;
extern int *jsort,*jsdsp,*jscnt,*jrdsp,*jrcnt;

extern void memoryuse();
extern void mpirange(int, int, int&, int&);

void mpiprej(int nj, int& mj) {
  int jsta,jend,i;

  jsort = new int [npmax];
  jsdsp = new int [nprocs];
  jscnt = new int [nprocs];
  jrdsp = new int [nprocs];
  jrcnt = new int [nprocs];
  mem = npmax*4+nprocs*4*4;
  memoryuse();

  mpirange(0,nj-1,jsta,jend);
  for( i=jsta; i<=jend; i++ ) {
    mj = i-jsta;
    xj[mj] = xj[i];
    yj[mj] = yj[i];
    zj[mj] = zj[i];
    gxj[mj] = gxj[i];
    gyj[mj] = gyj[i];
    gzj[mj] = gzj[i];
    vj[mj] = vj[i];
    sj[mj] = sj[i];
  }
  mj++;

}
