#include "mpi.h"
#include "../misc/constants.h"

extern float *xi,*yi,*zi,*gxi,*gyi,*gzi,*vi,*gxd,*gyd,*gzd,*vd;
extern int *nbi,**nxs,*nfn,*nek,*npart,*irank,*isort,*isdsp,*iscnt,*irdsp,*ircnt,*na,*nb;

extern void boxn(int, int, int);
extern void sort(int);
extern void sortvar(int, int, float*, int*);
extern void unsortvar(int, int, float*, int*);
extern void mpisortvar(int*, int*, int*, int*, float*);
extern void mpiunsortvar(int*, int*, int*, int*, float*);

void sorti(int& mi) {
  int lmaxd,nsubd,kl,i,j,iv,nsplit,lbi,nbc,ista,ic;
  double rb;

  lmaxd = lmax;
  while( int(pow(8,lmaxd)) < nprocs ) lmaxd++;
  nsubd = int(pow(8,lmaxd-lmax));
  kl = int(pow(2,lmaxd));
  rb = rd/kl;
  for( i=0; i<mi; i++ ) {
    nxs[0][i] = int((xi[i]-xmin)/rb);
    nxs[1][i] = int((yi[i]-ymin)/rb);
    nxs[2][i] = int((zi[i]-zmin)/rb);
    for( j=0; j<3; j++ ) {
      if( nxs[j][i] >= int(pow(2,lmaxd)) ) nxs[j][i] = nxs[j][i]-1;
    }
  }

  boxn(mi,3,lmaxd);
  for( i=0; i<mi; i++ ) {
    iv = nek[nfn[i]/nsubd]/nsub;
    nsplit = npart[iv+1]-npart[iv];
    na[i] = irank[npart[iv]]+nfn[i]%nsplit;
    nb[i] = i;
  }
  sort(mi);
  for( i=0; i<mi; i++ ) {
    isort[i] = nb[i];
  }

  lbi = 0;
  nbc = -1;
  for( i=0; i<nprocs; i++ ) iscnt[i] = 0;
  for( i=0; i<mi; i++ ) {
    if( na[i] != nbc ) {
      lbi++;
      if( lbi >= 2 ) {
        iscnt[na[i-1]] = i-ista;
      }
      ista = i;
      nbc = na[i];
    }
  }
  iscnt[na[mi-1]] = mi-ista;
  ic = 0;
  for( i=0; i<nprocs; i++ ) {
    isdsp[i] = ic;
    ic += iscnt[i];
  }

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Alltoall(iscnt,1,MPI_INT,ircnt,1,MPI_INT,MPI_COMM_WORLD);

  sortvar(0,mi,xi,isort);
  sortvar(0,mi,yi,isort);
  sortvar(0,mi,zi,isort);
  sortvar(0,mi,gxi,isort);
  sortvar(0,mi,gyi,isort);
  sortvar(0,mi,gzi,isort);
  sortvar(0,mi,vi,isort);

  mi = 0;
  for( i=0; i<nprocs; i++ ) {
    irdsp[i] = mi;
    mi += ircnt[i];
  }

  mpisortvar(isdsp,iscnt,irdsp,ircnt,xi);
  mpisortvar(isdsp,iscnt,irdsp,ircnt,yi);
  mpisortvar(isdsp,iscnt,irdsp,ircnt,zi);
  mpisortvar(isdsp,iscnt,irdsp,ircnt,gxi);
  mpisortvar(isdsp,iscnt,irdsp,ircnt,gyi);
  mpisortvar(isdsp,iscnt,irdsp,ircnt,gzi);
  mpisortvar(isdsp,iscnt,irdsp,ircnt,vi);

  for( i=0; i<mi; i++ ) {
    nxs[0][i] = int((xi[i]-xmin)/rb);
    nxs[1][i] = int((yi[i]-ymin)/rb);
    nxs[2][i] = int((zi[i]-zmin)/rb);
    for( j=0; j<3; j++ ) {
      if( nxs[j][i] >= int(pow(2,lmaxd)) ) nxs[j][i]--;
    }
  }

  boxn(mi,3,lmaxd);
  for( i=0; i<mi; i++ ) {
    na[i] = nfn[i];
    nb[i] = i;
  }
  sort(mi);
  for( i=0; i<mi; i++ ) nbi[i] = nb[i];

  sortvar(0,mi,xi,nbi);
  sortvar(0,mi,yi,nbi);
  sortvar(0,mi,zi,nbi);
  if( neq < 5 || 7 < neq ) {
    sortvar(0,mi,gxi,nbi);
    sortvar(0,mi,gyi,nbi);
    sortvar(0,mi,gzi,nbi);
  }
  sortvar(0,mi,vi,nbi);

}

void unsorti(int& mi) {
  int i;

  if( neq == 5 ) {
    for( i=0; i<npmax; i++ ) gxi[i] = 0;
    for( i=0; i<mi; i++ ) gxi[nbi[i]] = gxd[i];
  } else if( neq == 6 ) {
    for( i=0; i<npmax; i++ ) gyi[i] = 0;
    for( i=0; i<mi; i++ ) gyi[nbi[i]] = gxd[i];
  } else if( neq == 7 ) {
    for( i=0; i<npmax; i++ ) gzi[i] = 0;
    for( i=0; i<mi; i++ ) gzi[nbi[i]] = gxd[i];
  } else {
    for( i=0; i<npmax; i++ ) {
      gxi[i] = 0;
      gyi[i] = 0;
      gzi[i] = 0;
      vi[i] = 0;
    }
    for( i=0; i<mi; i++ ) {
      gxi[nbi[i]] = gxd[i];
      gyi[nbi[i]] = gyd[i];
      gzi[nbi[i]] = gzd[i];
      vi[nbi[i]] = vd[i];
    }
  }

  unsortvar(0,mi,xi,nbi);
  unsortvar(0,mi,yi,nbi);
  unsortvar(0,mi,zi,nbi);

  mpiunsortvar(isdsp,iscnt,irdsp,ircnt,xi);
  mpiunsortvar(isdsp,iscnt,irdsp,ircnt,yi);
  mpiunsortvar(isdsp,iscnt,irdsp,ircnt,zi);
  mpiunsortvar(isdsp,iscnt,irdsp,ircnt,gxi);
  mpiunsortvar(isdsp,iscnt,irdsp,ircnt,gyi);
  mpiunsortvar(isdsp,iscnt,irdsp,ircnt,gzi);
  mpiunsortvar(isdsp,iscnt,irdsp,ircnt,vi);

  mi = 0;
  for( i=0; i<nprocs; i++ ) mi += iscnt[i];

  unsortvar(0,mi,xi,isort);
  unsortvar(0,mi,yi,isort);
  unsortvar(0,mi,zi,isort);
  unsortvar(0,mi,gxi,isort);
  unsortvar(0,mi,gyi,isort);
  unsortvar(0,mi,gzi,isort);
  unsortvar(0,mi,vi,isort);

}
