#include "mpi.h"
#include "../misc/constants.h"

extern void memoryuse();
extern void memoryfree();

extern int *nfi,*nij,**neij,**nsij,**ncnt2,*nfrom,*nek,*npart,*irank,*ksdsp,*kscnt,*krdsp,*krcnt;

void jcnt(int& nr, int lbi, int lbj, int lev, int ipb) {
  int num,i,jj,ii,idig,nfdig,iv,inum,isrank,ij,ic,j;

  ncnt2 = new int* [nprocs];
  for( i=0; i<nprocs; i++ ) ncnt2[i] = new int [nbnp];
  mem = nprocs*nbnp*4;
  memoryuse();

  num = nprocs/pow(8,lev);
  if( num == 0 ) num = 1;
  for( i=0; i<nprocs; i++ ) kscnt[i] = 0;
  for( jj=0; jj<lbj; jj++ ) {
    for( i=0; i<nprocs; i++ ) {
      ncnt2[i][jj] = 0;
    }
  }
  for( ii=0; ii<lbi; ii++ ) {
    for( idig=0; idig<int(pow(8,lmax-lev)); idig++ ) {
      nfdig = nfi[ii]*int(pow(8,lmax-lev))+idig;
      if( nek[nfdig] != -1 ) {
        iv = nek[nfdig]/nsub;
        for( inum=npart[iv]; inum<npart[iv+1]; inum++ ) {
          isrank = irank[inum];
          if( isrank != myrank ) {
            for( ij=0; ij<nij[ii]; ij++ ) {
              jj = neij[ij][ii];
              if( ncnt2[isrank][jj] == 0 ) {
                nsij[kscnt[isrank]][isrank] = jj;
                ncnt2[isrank][jj] = 1;
                kscnt[isrank]++;
              }
            }
          }
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

  for( i=0; i<nprocs; i++ ) delete[] ncnt2[i];
  delete[] ncnt2;
  mem = nprocs*nbnp*4;
  memoryfree();
}
