#include "../misc/constants.h"

extern int *nfj,*nfl,*nfm,*nfo,*nc,*nd,**npx,*nnp;

extern void boxc(int, int, int*);
extern void jpsub(int, int, int, int, int&);

void jsbox(int lbj, int lbl, int lbm, int& lbjp, int lev, int ipb, int npb) {
  int ncmax,nmin,nmax,ntra,jc,jj,i,lb,k,ja,na,jb,nb,jx,jy,jz;
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
  }

  if( ntra == 1 ) {
    for( i=1; i<=3; i++ ) {
      if( i == 1 ) {
        lb = lbj;
        nd[2] = 0;
        for( jj=0; jj<lbj; jj++ ) nfo[jj] = nfj[jj];
      } else if ( i == 2 ) {
       lb = lbl;
       nd[2] = 1;
       for( jj=0; jj<lbl; jj++ ) nfo[jj] = nfl[jj];
      } else {
       lb = lbm;
       nd[2] = -1;
       for( jj=0; jj<lbm; jj++ ) nfo[jj] = nfm[jj];
      }
      for( jj=0; jj<lb; jj++ ) {
        boxc(nfo[jj],3,nc);
        for( k=0; k<2; k++ ) {
          if( nc[k] <= nmin ) {
            nd[k] = 1;
          } else if( nc[k] >= nmax ) {
            nd[k] = -1;
          } else {
            nd[k] = 0;
          }
        }
        if( i == 1 || (i == 2 && nc[2] <= nmin) || (i == 3 && nc[2] >= nmax) ) {
          if( nd[2] != 0 ) jpsub(0,0,nd[2],jj,jc);
          if( nd[1] != 0 ) jpsub(0,nd[1],nd[2],jj,jc);
          if( nd[0] != 0 ) jpsub(nd[0],0,nd[2],jj,jc);
          if( nd[0]*nd[1] != 0 ) jpsub(nd[0],nd[1],nd[2],jj,jc);
        }
      }
    }
  } else if( ntra == 2 ) {
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
      for( ja=0; ja<3; ja++ ) {
        na = ja-1;
        for( jb=0; jb<3; jb++ ) {
          nb = jb-1;
          if( nd[0] != 0 ) jpsub(nd[0],na,nb,jj,jc);
          if( nd[1] != 0 ) jpsub(na,nd[1],nb,jj,jc);
          if( nd[2] != 0 ) jpsub(na,nb,nd[2],jj,jc);
        }
        if( nd[0]*nd[1] != 0 ) jpsub(nd[0],nd[1],na,jj,jc);
        if( nd[1]*nd[2] != 0 ) jpsub(na,nd[1],nd[2],jj,jc);
        if( nd[2]*nd[0] != 0 ) jpsub(nd[0],na,nd[2],jj,jc);
      }
      if( nd[0]*nd[1]*nd[2] != 0 ) jpsub(nd[0],nd[1],nd[2],jj,jc);
      for( jx=0; jx<3; jx++ ) {
        nd[0] = jx-1;
        for( jy=0; jy<3; jy++ ) {
          nd[1] = jy-1;
          for( jz=0; jz<3; jz++ ) {
            nd[2] = jz-1;
            if (jx*jy*jz != 1) jpsub(nd[0],nd[1],nd[2],jj,jc);
          }
        }
      }
    }
  }
  lbjp = jc;
}
