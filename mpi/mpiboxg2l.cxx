#include "mpi.h"
#include "../misc/constants.h"

extern int *nxadj,*nadjncy,*nvwgt,*nadjwgt,*nvtxdist,*mxadj,*madjncy,*mvwgt,*madjwgt;
extern int *ksdsp,*kscnt,*krdsp,*krcnt,*lsdsp,*lscnt,*lrdsp,*lrcnt;

extern void memoryuse();
extern void memoryfree();

void boxg2l(int& nv) {
  int i,iwork1,iwork2,ista,iend,ic,nadj;

  mxadj = new int [nbne+1];
  madjncy = new int [npmax];
  mvwgt = new int [nbne];
  madjwgt = new int [npmax];
  ksdsp = new int [nprocs];
  kscnt = new int [nprocs];
  lsdsp = new int [nprocs];
  lscnt = new int [nprocs];
  lrdsp = new int [nprocs];
  lrcnt = new int [nprocs];
  mem = nbne*2*4+npmax*2*4+nprocs*6*4;
  memoryuse();

  for( i=0; i<nprocs; i++ ) {
    iwork1 = nv/nprocs;
    iwork2 = nv%nprocs;
    ista = i*iwork1+std::min(i,iwork2);
    iend = ista+iwork1;
    if( iwork2 > i ) iend++;
    ksdsp[i] = ista;
    kscnt[i] = iend-ista;
    lsdsp[i] = nxadj[ista];
    lscnt[i] = nxadj[iend];
  }
  if( nv > nprocs ) {
    MPI_Alltoall(kscnt,1,MPI_INT,krcnt,1,MPI_INT,MPI_COMM_WORLD);
    MPI_Alltoall(lscnt,1,MPI_INT,lrcnt,1,MPI_INT,MPI_COMM_WORLD);
  } else {
    for( i=0; i<nprocs; i++ ) {
      krcnt[i] = kscnt[i];
      lrcnt[i] = lscnt[i];
    }
  }

  ic = 0;
  nvtxdist[0] = 0;
  for( i=0; i<nprocs; i++ ) {
    krdsp[i] = ic;
    ic += krcnt[i];
    nvtxdist[i+1] = ic;
  }
  MPI_Scatter(lsdsp,1,MPI_INT,&ista,1,MPI_INT,0,MPI_COMM_WORLD);
  nv = krcnt[myrank];
  MPI_Scatterv(nxadj+1,kscnt,ksdsp,MPI_INT,mxadj+1,nv,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Scatterv(nvwgt,kscnt,ksdsp,MPI_INT,mvwgt,nv,MPI_INT,0,MPI_COMM_WORLD);
  nadj = lrcnt[myrank];
  MPI_Scatterv(nadjncy,lscnt,lsdsp,MPI_INT,madjncy,nadj,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Scatterv(nadjwgt,lscnt,lsdsp,MPI_INT,madjwgt,nadj,MPI_INT,0,MPI_COMM_WORLD);
  nxadj[0] = 0;
  for( i=0; i<nv; i++ ) {
    nxadj[i+1] = mxadj[i+1]-ista;
    nvwgt[i] = mvwgt[i];
  }
  for( i=0; i<nadj; i++ ) {
    nadjncy[i] = madjncy[i];
    nadjwgt[i] = madjwgt[i];
  }

  delete[] mxadj;
  delete[] madjncy;
  delete[] mvwgt;
  delete[] madjwgt;
  delete[] ksdsp;
  delete[] kscnt;
  delete[] lsdsp;
  delete[] lscnt;
  delete[] lrdsp;
  delete[] lrcnt;
  mem = nbne*2*4+npmax*2*4+nprocs*6*4;
  memoryfree();
}
