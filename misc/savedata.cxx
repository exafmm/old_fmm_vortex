#include "../misc/constants.h"

extern float *xj,*yj,*zj,*gxj,*gyj,*gzj,*vj,*sj;
extern double **ptr;
extern int **ntr,*nne;

void savedata(int nwe, int np, int nge) {
  int nnd,nbyte,i,j,ij,nwn;
  int nij[6*nwmax];
  std::fstream fid;
  fid.open("particles.dat",std::ios::out | std::ios::binary | std::ios::ate);

  nnd = 6;
  
  nbyte = sizeof(int);
  fid.write((char *)(&nbyte),sizeof(int));
  fid.write((char *)(&np),sizeof(int));
  fid.write((char *)(&nbyte),sizeof(int));
  if( nwe != 0 && nge != 8 ) {
    for( i=0; i<nwe; i++ ) {
      for( j=0; j<nnd; j++ ) {
        ij = i*nnd+j;
        nij[ij] = ntr[j][i];
      }
    }
    nwn = nwe*nnd;
    nbyte = sizeof(double)*nwn*3;
    fid.write((char *)(&nbyte),sizeof(int));
    for( i=0; i<nwn; i++ ) fid.write((char *)(&ptr[0][nij[i]]),sizeof(double));
    for( i=0; i<nwn; i++ ) fid.write((char *)(&ptr[1][nij[i]]),sizeof(double));
    for( i=0; i<nwn; i++ ) fid.write((char *)(&ptr[2][nij[i]]),sizeof(double));
    fid.write((char *)(&nbyte),sizeof(int));
  }

  nbyte = sizeof(float)*np*8;
  fid.write((char *)(&nbyte),sizeof(int));
  for( i=0; i<np; i++ ) fid.write((char *)(&xj[i]),sizeof(float));
  for( i=0; i<np; i++ ) fid.write((char *)(&yj[i]),sizeof(float));
  for( i=0; i<np; i++ ) fid.write((char *)(&zj[i]),sizeof(float));
  for( i=0; i<np; i++ ) fid.write((char *)(&gxj[i]),sizeof(float));
  for( i=0; i<np; i++ ) fid.write((char *)(&gyj[i]),sizeof(float));
  for( i=0; i<np; i++ ) fid.write((char *)(&gzj[i]),sizeof(float));
  for( i=0; i<np; i++ ) fid.write((char *)(&vj[i]),sizeof(float));
  for( i=0; i<np; i++ ) fid.write((char *)(&sj[i]),sizeof(float));
  fid.write((char *)(&nbyte),sizeof(int));
  fid.close();
}
