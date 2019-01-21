#include "mpi.h"
#include "parmetis.h"
#include "../misc/constants.h"

extern int *nfi,*nfj,**npx,*nek,*ixadj,*nxadj,*nadjncy,*nvwgt,*nadjwgt,*nvtxdist;
extern int *npart,*npartd,*irank,*krdsp,*krcnt;

extern void memoryuse();
extern void memoryfree();
extern void boxdatai(int, int, int, int&, double&);
extern void boxdataj(int, int, int, int&, double&);
extern void boxl2g(int&, int&);
extern void ijbox(int, int, int, int, int);
extern void setvert(int, int);
extern void setedge(int, int);
extern void boxg2l(int&);

void boxpart(int mi, int mj, int lev) {
  int lbi,lbj,lbk,ii,jj,nv,iv,i,j,nvd,iwork1,iwork2,ista,iend;
  int nwgtflag=3,numflag=0,ncon=1,nparts,noptions[5]={0,0,0,0,0},nedgecut;
  float tpwgts[nprocs],ubvec=1.05;
  double rb;
  MPI_Comm comm = MPI_COMM_WORLD;

  krdsp = new int [nprocs];
  krcnt = new int [nprocs];
  mem = nprocs*2*4;
  memoryuse();

  nsub = int(pow(8,lmax))/nprocs/64;
  if( nsub < 1 ) nsub = 1;
  boxdatai(1,mi,lev,lbi,rb);
  boxdataj(1,mj,lev,lbj,rb);

#if 0
  boxl2g(lbi,lbj);
  lbk = 0;
  for( ii=0; ii<nbmax; ii++ ) nek[ii] = -1;
  for( ii=0; ii<lbi; ii++ ) {
    nek[nfi[ii]] = lbk;
    lbk++;
  }
  for( jj=0; jj<lbj; jj++ ) {
    if( nek[nfj[jj]] == -1 ) {
      nek[nfj[jj]] = lbk;
      lbk++;
    }
  }
  nv = (lbk-1)/nsub+1;
  nparts = std::min(nv,nprocs);
  for( iv=0; iv<nv; iv++ ) {
    nxadj[iv] = 0;
    nvwgt[iv] = 0;
  }
  nxadj[nv] = 0;
  for( i=0; i<npmax; i++ ) {
    nadjncy[i] = 0;
    nadjwgt[i] = 0;
  }
  for( i=0; i<nbnp; i++ ) {
    for( j=0; j<3; j++ ) {
      npx[j][i] = 0;
    }
  }
  ijbox(lbi,lbj,lev,-2,0);
  setvert(lbi,2);
  ijbox(lbi,lbj,lev,-1,0);
  setvert(lbi,4);
  for( iv=0; iv<nv; iv++ ) {
    ixadj[iv] = nxadj[iv];
    nxadj[iv+1] += nxadj[iv];
  }
  ijbox(lbi,lbj,lev,-2,0);
  setedge(lbi,2);
  ijbox(lbi,lbj,lev,-1,0);
  setedge(lbi,4);
  boxg2l(nv);

  for( i=0; i<nparts; i++ ) tpwgts[i] = 1.0/nparts;
  ParMETIS_V3_PartKway(nvtxdist,nxadj,nadjncy,nvwgt,nadjwgt,&nwgtflag,&numflag,&ncon,&nparts,
                       tpwgts,&ubvec,noptions,&nedgecut,npart,&comm);
  MPI_Allgather(&nv,1,MPI_INT,krcnt,1,MPI_INT,MPI_COMM_WORLD);
  nvd = 0;
  for( i=0; i<nprocs; i++ ) {
    krdsp[i] = nvd;
    nvd += krcnt[i];
  }

  MPI_Allgatherv(npart,nv,MPI_INT,npartd,krcnt,krdsp,MPI_INT,MPI_COMM_WORLD);
#else
  boxl2g(lbi,lbj);
  lbk = 0;
  for( ii=0; ii<nbmax; ii++ ) nek[ii] = -1;
  for( ii=0; ii<lbi; ii++ ) {
    nek[nfi[ii]] = lbk;
    lbk++;
  }
  for( jj=0; jj<lbj; jj++ ) {
    if( nek[nfj[jj]] == -1 ) {
      nek[nfj[jj]] = lbk;
      lbk++;
    }
  }
  nvd = (lbk-1)/nsub+1;
  nparts = std::min(nvd,nprocs);
#endif

  npart[0] = 0;
  nv = nvd;
  for( i=0; i<nv; i++ ) {
    iwork1 = nprocs/nv;
    iwork2 = nprocs%nv;
    if( nprocs < nv ) iwork2 = nv;
    ista = i*iwork1+std::min(i,iwork2);
    iend = ista+iwork1;
    if( iwork2 > i ) iend++;
    for( j=ista; j<iend; j++ ) {
      irank[j] = j/(nv/nparts); // Morton order partition
//      irank[j] = j%nprocs; // worst partition
//      irank[j] = npartd[i]; // parmetis partition
    }
    npart[i+1] = iend;
  }

  delete[] krdsp;
  delete[] krcnt;
  mem = nprocs*2*4;
  memoryfree();
}
