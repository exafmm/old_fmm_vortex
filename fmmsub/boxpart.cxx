#include "../misc/constants.h"

extern int *nfi,*nfj,*nek,*npart,*irank;

extern void boxdatai(int, int, int, int&, double&);
extern void boxdataj(int, int, int, int&, double&);

void boxpart(int mi, int mj, int lev) {
  int lbi,lbj,lbk,ii,jj,nv,i,iwork1,iwork2,ista,iend,j;
  double rb;

  boxdatai(1,mi,lev,lbi,rb);
  boxdataj(1,mj,lev,lbj,rb);
  lbk = 0;
  for( ii=0; ii<nbmax; ii++ ) nek[ii] = 0;
  for( ii=0; ii<lbi; ii++ ) {
    nek[nfi[ii]] = lbk;
    lbk++;
  }
  for( jj=0; jj<lbj; jj++ ) {
    if( nek[nfj[jj]] == 0 ) {
      nek[nfj[jj]] = lbk;
      lbk++;
    }
  }
  nv = lbk;

  npart[0] = 0;
  for( i=0; i<nv; i++ ) {
    iwork1 = nprocs/nv;
    iwork2 = nprocs%nv;
    if( nprocs < nv ) iwork2 = nv;
    ista = i*iwork1+std::min(i,iwork2);
    iend = ista+iwork1;
    if( iwork2 > i ) iend++;
    for( j=ista; j<iend; j++ ) irank[j] = 0;
    npart[i+1] = iend;
  }

}
