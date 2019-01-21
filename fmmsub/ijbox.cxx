#include "../misc/constants.h"

extern int *nei,*nfi,*nej,*nfj,*nlbj,*nc,**npx,**neij,*nij,*njb;

extern void boxc(int, int, int*);
extern void boxn1(int*, int&, int);

void ijbox(int lbi, int lbj, int lev, int ipb, int npb) {
  int nmin,neib,ixmin,ixmax,iymin,iymax,izmin,izmax,ii,jj,jx,jy,jz,ix,iy,iz,ie;
  int jxp,jyp,jzp,ixp,iyp,izp;

  if( ipb == -3 ) {
    nmin = 0;
    neib = 3;
  } else if( ipb == -2 ) {
    nmin = 1;
    neib = 2;
  } else if( ipb == -1 ) {
    nmin = 2;
    neib = 4;
  } else if( ipb == npb ) {
    nmin = 1;
    neib = 3;
  } else {
    nmin = 3;
    neib = 4;
  }
  ixmin = 1000000;
  ixmax = -1000000;
  iymin = 1000000;
  iymax = -1000000;
  izmin = 1000000;
  izmax = -1000000;
  for( ii=0; ii<lbi; ii++ ) {
    nij[ii] = 0;
    boxc(nfi[ii],3,nc);
    ixmin = std::min(ixmin,nc[0]);
    ixmax = std::max(ixmax,nc[0]);
    iymin = std::min(iymin,nc[1]);
    iymax = std::max(iymax,nc[1]);
    izmin = std::min(izmin,nc[2]);
    izmax = std::max(izmax,nc[2]);
  }
  if( neib == 2 ) {
    for( jj=0; jj<lbj; jj++ ) {
      boxc(nfj[jj],3,nc);
      jx = nc[0]+npx[0][jj]*int(pow(2,lev));
      jy = nc[1]+npx[1][jj]*int(pow(2,lev));
      jz = nc[2]+npx[2][jj]*int(pow(2,lev));
      for( ix=std::max(jx-1,ixmin); ix<=std::min(jx+1,ixmax); ix++ ) {
        for( iy=std::max(jy-1,iymin); iy<=std::min(jy+1,iymax); iy++ ) {
          for( iz=std::max(jz-1,izmin); iz<=std::min(jz+1,izmax); iz++ ) {
            nc[0] = ix;
            nc[1] = iy;
            nc[2] = iz;
            boxn1(nc,ie,lev);
            ii = nei[ie];
            if( ii != -1 ) {
              neij[nij[ii]][ii] = jj;
              nij[ii]++;
            }
          }
        }
      }
    }
  } else if( neib == 3 ) {
    for( jj=0; jj<lbj; jj++ ) {
      boxc(nfj[jj],3,nc);
      jx = nc[0]+npx[0][jj]*int(pow(2,lev));
      jy = nc[1]+npx[1][jj]*int(pow(2,lev));
      jz = nc[2]+npx[2][jj]*int(pow(2,lev));
      jxp = (jx+nmin)/2;
      jyp = (jy+nmin)/2;
      jzp = (jz+nmin)/2;
      for( ii=0; ii<lbi; ii++ ) {
        boxc(nfi[ii],3,nc);
        ix = nc[0];
        iy = nc[1];
        iz = nc[2];
        ixp = (ix+nmin)/2;
        iyp = (iy+nmin)/2;
        izp = (iz+nmin)/2;
        if( ix < jx-1 || jx+1 < ix || iy < jy-1 || jy+1 < iy || iz < jz-1 || jz+1 < iz ) {
          neij[nij[ii]][ii] = jj;
          nij[ii]++;
        }
      }
    }
  } else if( neib == 4 ) {
    for( jj=0; jj<lbj; jj++ ) {
      boxc(nfj[jj],3,nc);
      jx = nc[0]+npx[0][jj]*int(pow(2,lev));
      jy = nc[1]+npx[1][jj]*int(pow(2,lev));
      jz = nc[2]+npx[2][jj]*int(pow(2,lev));
      jxp = (jx+nmin)/2;
      jyp = (jy+nmin)/2;
      jzp = (jz+nmin)/2;
      for( ixp=jxp-1; ixp<=jxp+1; ixp++ ) {
        for( iyp=jyp-1; iyp<=jyp+1; iyp++ ) {
          for( izp=jzp-1; izp<=jzp+1; izp++ ) {
            for( ix=std::max(2*ixp-nmin,ixmin); ix<=std::min(2*ixp-nmin+1,ixmax); ix++ ) {
              for( iy=std::max(2*iyp-nmin,iymin); iy<=std::min(2*iyp-nmin+1,iymax); iy++ ) {
                for( iz=std::max(2*izp-nmin,izmin); iz<=std::min(2*izp-nmin+1,izmax); iz++ ) {
                  if( ix < jx-1 || jx+1 < ix || iy < jy-1 || jy+1 < iy || iz < jz-1 || jz+1 < iz ) {
                    nc[0] = ix;
                    nc[1] = iy;
                    nc[2] = iz;
                    boxn1(nc,ie,lev);
                    ii = nei[ie];
                    if( ii != -1 ) {
                      neij[nij[ii]][ii] = jj;
                      nij[ii]++;
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  for( jj=0; jj<lbj; jj++ ) {
    njb[jj] = nej[nfj[jj]]+nlbj[lev-1];
  }
  if( ipb > 0 ) {
    for( jj=0; jj<lbj; jj++ ) {
      njb[jj] = ipb+nbnet-1;
    }
  }
}
