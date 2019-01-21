#include "mpi.h"
#include "../misc/constants.h"

extern int *nfi,*nij,**neij,**nsij,*nform,*nek,*npart,*irank,*ksdsp,*kscnt,*krdsp,*krcnt;

void jcnt(int& nr, int lbi, int lbj, int lev, int ipb) {
  int num,nom,i,jj,ij,iv,inum,isrank,ic,ncnt[nprocs];

  num = nprocs/int(pow(8,lev));
  if( num == 0 ) num = 1;
  for( i=0; i<nprocs; i++ ) kscnt[i] = 0;
  for( jj=0; jj<lbj; jj++ ) {
    for( i=0; i<nprocs; i++ ) ncnt[i] = 0;
    for( ij=0; ij<nij[jj]; ij++ ) {
      if( nek[nfi[ii]] != -1 ) {
        iv = nek[nfi[i]]/nsub;
        for( inum=npart[iv]; inum<npart[iv+1]; inum++ ) {
        isrank = irank[inum];
        if( isrank != myrank && ncnt[isrank] == 0 ) {
          nsij[kscnt[isrank]][isrank] = jj;
          ncnt[isrank] = 1;
          kscnt[isrank]++;
        }
      }
    }
  }
  if( myrank%num != 0 && ( ipb == -1 || ipb == -3 ) ) {
    for( i=0; i<nprocs; i++ ) kscnt[i] = 0;
  }

  ic = 0;
  for( i=0; i<nprocs; i++ ) {
    ksdsp[i] = ic;
    ic += kscnt[i];
  }

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Alltoall(kscnt,1,MPI_INT,krcnt,1,MPI_INT,MPI_COMM_WORLD);

  nr = 0;
  for( i=0; i<nprocs; i++ ) {
    krdsp[i] = nr;
    nr += krcnt[i];
    for( j=krdsp[i]; j<krdsp[i]+krcnt[i]; j++ ) nfrom[j] = i;
  }

}
