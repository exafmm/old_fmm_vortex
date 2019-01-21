#include "../misc/constants.h"

extern int *nfj,*nfo,*nc,*nd,**npx,*nnp;

extern void boxc(int, int, int*);
extern void jpsub(int, int, int, int, int&);

void jcbox(int lbj, int& lbjp, int lev, int ipb, int npb) {
  int ncmax,nmin,nmax,ntra,jc,jj,i,k,ja,na,jx,jy;

  ncmax = int(pow(2,lev))-1;
  if( ipb == -1 ) {
    nmin = 1;
    nmax = ncmax-1;
    ntra = 1;
  } else if( ipb == npb || ipb == -2 ) {
    nmin = 0;
    nmax = ncmax;
    ntra = 1;
  } else {
    nmin = 0;
    nmax = ncmax;
    ntra = 2;
  }
  jc = lbj;
  for( jj=0; jj<lbj; jj++ ) {
    for( i=0; i<3; i++ ) npx[i][jj] = 0;
    nnp[jj] = jj;
    nfo[jj] = nfj[jj];
  }
  for( jj=0; jj<lbj; jj++ ) {
    boxc(nfj[jj],3,nc);
    for( k=0; k<3; k++ ) {
      if( nc[k] <= nmin ) {
        nd[k] = ntra;
      } else if( nc[k] >= nmax ) {
        nd[k] = -ntra;
      } else {
        nd[k] = 0;
      }
    }
    if( ntra == 1 ) {
      if( nd[0] != 0 ) jpsub(nd[0],0,0,jj,jc);
      if( nd[1] != 0 ) jpsub(0,nd[1],0,jj,jc);
      if( nd[0]*nd[1] != 0 ) jpsub(nd[0],nd[1],0,jj,jc);
    } else if ( ntra == 2 ) {
      for( ja=0; ja<3; ja++ ) {
        na = ja-1;
        if( nd[0] != 0 ) jpsub(nd[0],na,0,jj,jc);
        if( nd[1] != 0 ) jpsub(na,nd[1],0,jj,jc);
      }
      if( nd[0]*nd[1] != 0 ) jpsub(nd[0],nd[1],0,jj,jc);
      for( jx=0; jx<3; jx++ ) {
        nd[0] = jx-1;
        for( jy=0; jy<3; jy++ ) {
          nd[1] = jy-1;
          if (jx*jy != 1) jpsub(nd[0],nd[1],0,jj,jc);
        }
      }
    }
  }
  lbjp = jc;
}
