#include "mpi.h"
#include "../misc/parameters.h"
#include "../misc/constants.h"

void nlevel(int mmax){
  int nmax,ned,mprocs,mpi,ndf,ienv;
  int nlev[1000000];
  std::ifstream fid;

  MPI_Allreduce(&mmax,&nmax,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  ned = 0;
  if( neq < 0 ) {
    ned = neq;
    neq = 0;
  }
  fid.open("../../dat/nlevel",std::ios::binary);
  mprocs = nprocs;
  mpi = -1;
  while ( mprocs != 0 ) {
    mpi++;
    mprocs /= 2;
  }
  ndf = ndev*2+nfmm%2;

  if( neq < 0 ) {
    fid.read((char *)(&nlev),4);
    fid.read((char *)(&nlev),sizeof(nlev));
    ienv = 100000*mpi+10000*ndf+1000*nge+100*nsmt-10*neq;
    nlev[ienv-neq-2] = nmax;
    fid.seekg(0,std::ios::beg);
  } else if( lmax >= 0 ) {
    fid.read((char *)(&nlev),4);
    fid.read((char *)(&nlev),sizeof(nlev));
    ienv = 100000*mpi+10000*ndf+1000*nge+100*nsmt+10*neq;
    lmax = 1;
    if( nmax < nlev[ienv] ) {
      lmax += 1;
    } else if( nmax < nlev[ienv+1] ) {
      lmax += 2;
    } else if( nmax < nlev[ienv+2] ) {
      lmax += 3;
    } else if( nmax < nlev[ienv+3] ) {
      lmax += 4;
    } else if( nmax < nlev[ienv+4] ) {
      lmax += 5;
    } else if( nmax < nlev[ienv+5] ) {
      lmax += 6;
    } else {
      lmax += 7;
    }
  } else {
    lmax = -lmax;
  }
  if( ned != 0 ) neq = ned;
  fid.close();
}
