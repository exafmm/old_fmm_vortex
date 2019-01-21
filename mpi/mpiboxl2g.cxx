#include "mpi.h"
#include "../misc/constants.h"

extern int **ndi,*nei,*nfi,**ndj,*nej,*nfj,**ndo,*krdsp,*krcnt,*na,*nb;

extern void memoryuse();
extern void memoryfree();
extern void sort(int);

void boxl2g(int& lbi, int& lbj) {
  int i,lbid,nbc,lbjd;

  ndo = new int* [2];
  for( i=0; i<2; i++ ) ndo[i] = new int [npmax];
  mem = npmax*2*4;
  memoryuse();

  for( i=0; i<nbmax; i++ ) {
    nei[i] = 0;
    nej[i] = 0;
  }

  MPI_Allgather(&lbi,1,MPI_INT,krcnt,1,MPI_INT,MPI_COMM_WORLD);
  lbid = 0;
  for( i=0; i<nprocs; i++ ) {
    krdsp[i] = lbid;
    lbid += krcnt[i];
  }

  MPI_Allgatherv(nfi,lbi,MPI_INT,na,krcnt,krdsp,MPI_INT,MPI_COMM_WORLD);
  MPI_Allgatherv(ndi[0],lbi,MPI_INT,ndo[0],krcnt,krdsp,MPI_INT,MPI_COMM_WORLD);
  MPI_Allgatherv(ndi[1],lbi,MPI_INT,ndo[1],krcnt,krdsp,MPI_INT,MPI_COMM_WORLD);
  sort(lbid);

  lbi = 0;
  nbc = -1;
  for( i=0; i<lbid; i++ ) {
    if( na[i] != nbc ) {
      nei[na[i]] = lbi;
      nfi[lbi] = na[i];
      ndi[0][lbi] = ndo[0][i];
      ndi[1][lbi] = ndo[1][i];
      nbc = na[i];
      lbi++;
    }
  }

  MPI_Allgather(&lbj,1,MPI_INT,krcnt,1,MPI_INT,MPI_COMM_WORLD);
  lbjd = 0;
  for( i=0; i<nprocs; i++ ) {
    krdsp[i] = lbjd;
    lbjd += krcnt[i];
  }
  MPI_Allgatherv(nfj,lbj,MPI_INT,na,krcnt,krdsp,MPI_INT,MPI_COMM_WORLD);
  MPI_Allgatherv(ndj[0],lbj,MPI_INT,ndo[0],krcnt,krdsp,MPI_INT,MPI_COMM_WORLD);
  MPI_Allgatherv(ndj[1],lbj,MPI_INT,ndo[1],krcnt,krdsp,MPI_INT,MPI_COMM_WORLD);
  sort(lbjd);
  lbj = 0;
  nbc = -1;
  for( i=0; i<lbjd; i++ ) {
    if( na[i] != nbc ) {
      nej[na[i]] = lbj;
      nfj[lbj] = na[i];
      ndj[0][lbj] = ndo[0][i];
      ndj[1][lbj] = ndo[1][i];
      nbc = na[i];
      lbj++;
    }
  }

  for( i=0; i<2; i++ ) delete[] ndo[i];
  delete[] ndo;
  mem = npmax*2*4;
  memoryfree();
}
