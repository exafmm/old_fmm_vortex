#include "mpi.h"
#include "../misc/constants.h"

extern float *xj,*yj,*zj,*gxj,*gyj,*gzj,*vj,*sj;
extern int *nbj,**nxs,*nfn,*nek,*npart,*irank,*jsort,*jsdsp,*jscnt,*jrdsp,*jrcnt,*na,*nb;

extern void boxn(int, int, int);
extern void sort(int);
extern void sortvar(int, int, float*, int*);
extern void unsortvar(int, int, float*, int*);
extern void mpisortvar(int*, int*, int*, int*, float*);
extern void mpiunsortvar(int*, int*, int*, int*, float*);

void sortj(int& mj) {
  int lmaxd,nsubd,kl,i,j,jv,nsplit,lbj,nbc,ista,ic;
  double rb;

  lmaxd = lmax;
  while( int(pow(8,lmaxd)) < nprocs ) lmaxd++;
  nsubd = int(pow(8,lmaxd-lmax));
  kl = int(pow(2,lmaxd));
  rb = rd/kl;
  for( i=0; i<mj; i++ ) {
    nxs[0][i] = int((xj[i]-xmin)/rb);
    nxs[1][i] = int((yj[i]-ymin)/rb);
    nxs[2][i] = int((zj[i]-zmin)/rb);
    for( j=0; j<3; j++ ) {
      if( nxs[j][i] >= int(pow(2,lmaxd)) ) nxs[j][i] = nxs[j][i]-1;
    }
  }

  boxn(mj,3,lmaxd);
  for( i=0; i<mj; i++ ) {
    jv = nek[nfn[i]/nsubd]/nsub;
    nsplit = npart[jv+1]-npart[jv];
    na[i] = irank[npart[jv]]+nfn[i]%nsplit;
    nb[i] = i;
  }
  sort(mj);
  for( i=0; i<mj; i++ ) {
    jsort[i] = nb[i];
  }

  lbj = 0;
  nbc = -1;
  for( i=0; i<nprocs; i++ ) jscnt[i] = 0;
  for( i=0; i<mj; i++ ) {
    if( na[i] != nbc ) {
      lbj++;
      if( lbj >= 2 ) {
        jscnt[na[i-1]] = i-ista;
      }
      ista = i;
      nbc = na[i];
    }
  }
  jscnt[na[mj-1]] = mj-ista;
  ic = 0;
  for( i=0; i<nprocs; i++ ) {
    jsdsp[i] = ic;
    ic += jscnt[i];
  }

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Alltoall(jscnt,1,MPI_INT,jrcnt,1,MPI_INT,MPI_COMM_WORLD);

  sortvar(0,mj,xj,jsort);
  sortvar(0,mj,yj,jsort);
  sortvar(0,mj,zj,jsort);
  sortvar(0,mj,gxj,jsort);
  sortvar(0,mj,gyj,jsort);
  sortvar(0,mj,gzj,jsort);
  sortvar(0,mj,vj,jsort);
  sortvar(0,mj,sj,jsort);

  mj = 0;
  for( i=0; i<nprocs; i++ ) {
    jrdsp[i] = mj;
    mj += jrcnt[i];
  }

  mpisortvar(jsdsp,jscnt,jrdsp,jrcnt,xj);
  mpisortvar(jsdsp,jscnt,jrdsp,jrcnt,yj);
  mpisortvar(jsdsp,jscnt,jrdsp,jrcnt,zj);
  mpisortvar(jsdsp,jscnt,jrdsp,jrcnt,gxj);
  mpisortvar(jsdsp,jscnt,jrdsp,jrcnt,gyj);
  mpisortvar(jsdsp,jscnt,jrdsp,jrcnt,gzj);
  mpisortvar(jsdsp,jscnt,jrdsp,jrcnt,vj);
  mpisortvar(jsdsp,jscnt,jrdsp,jrcnt,sj);

  for( i=0; i<mj; i++ ) {
    nxs[0][i] = int((xj[i]-xmin)/rb);
    nxs[1][i] = int((yj[i]-ymin)/rb);
    nxs[2][i] = int((zj[i]-zmin)/rb);
    for( j=0; j<3; j++ ) {
      if( nxs[j][i] >= int(pow(2,lmaxd)) ) nxs[j][i]--;
    }
  }

  boxn(mj,3,lmaxd);
  for( i=0; i<mj; i++ ) {
    na[i] = nfn[i];
    nb[i] = i;
  }
  sort(mj);
  for( i=0; i<mj; i++ ) nbj[i] = nb[i];

  sortvar(0,mj,xj,nbj);
  sortvar(0,mj,yj,nbj);
  sortvar(0,mj,zj,nbj);
  sortvar(0,mj,gxj,nbj);
  sortvar(0,mj,gyj,nbj);
  sortvar(0,mj,gzj,nbj);
  sortvar(0,mj,vj,nbj);
  sortvar(0,mj,sj,nbj);

}

void unsortj(int& mj) {
  int i;

  unsortvar(0,mj,xj,nbj);
  unsortvar(0,mj,yj,nbj);
  unsortvar(0,mj,zj,nbj);
  unsortvar(0,mj,gxj,nbj);
  unsortvar(0,mj,gyj,nbj);
  unsortvar(0,mj,gzj,nbj);
  unsortvar(0,mj,vj,nbj);
  unsortvar(0,mj,sj,nbj);

  mpiunsortvar(jsdsp,jscnt,jrdsp,jrcnt,xj);
  mpiunsortvar(jsdsp,jscnt,jrdsp,jrcnt,yj);
  mpiunsortvar(jsdsp,jscnt,jrdsp,jrcnt,zj);
  mpiunsortvar(jsdsp,jscnt,jrdsp,jrcnt,gxj);
  mpiunsortvar(jsdsp,jscnt,jrdsp,jrcnt,gyj);
  mpiunsortvar(jsdsp,jscnt,jrdsp,jrcnt,gzj);
  mpiunsortvar(jsdsp,jscnt,jrdsp,jrcnt,vj);
  mpiunsortvar(jsdsp,jscnt,jrdsp,jrcnt,sj);

  mj = 0;
  for( i=0; i<nprocs; i++ ) mj += jscnt[i];

  unsortvar(0,mj,xj,jsort);
  unsortvar(0,mj,yj,jsort);
  unsortvar(0,mj,zj,jsort);
  unsortvar(0,mj,gxj,jsort);
  unsortvar(0,mj,gyj,jsort);
  unsortvar(0,mj,gzj,jsort);
  unsortvar(0,mj,vj,jsort);
  unsortvar(0,mj,sj,jsort);

}
