#include "../misc/constants.h"

extern int *nc;
extern double (*gx)[nspm],(*gy)[nspm],(*gz)[nspm];

extern void boxc(int, int, int*);
extern void pm2m(double*, double*, double*, double*, double*, double*, int, int);

void m2m2(int nmp, int mp, int ip, int jp, int nc3, double rb) {
  int mc3,i,j,ib,jb,nfjc;
  double ggx[nspm],ggy[nspm],ggz[nspm],ggxd[nspm],ggyd[nspm],ggzd[nspm];

  mc3 = nc3;
  if( 2 <= nc3 && nc3 <= 7 ) {
    for( i=0; i<nmp; i++ ) {
      gx[ip-1][i] = 0;
      gy[ip-1][i] = 0;
      gz[ip-1][i] = 0;
    }
    mc3 = nc3-2;
  } else if ( nc3 == 9 ) {
    for( i=0; i<8; i++ ) {
      for( j=0; j<nmp; j++ ) {
        gx[i+3*nbnet][j] = 0;
        gy[i+3*nbnet][j] = 0;
        gz[i+3*nbnet][j] = 0;
      }
    }
  }
  for( i=0; i<8; i++ ) {
    boxc(i,3,nc);
// M2M @ lev 1 to sublev 1
    if( nc[2] == mc3 ) {
      ib = ip-1;
      jb = i+jp;
      nfjc = 7-i;
      for( j=0; j<nmp; j++ ) {
        ggxd[j] = gx[jb][j];
        ggyd[j] = gy[jb][j];
        ggzd[j] = gz[jb][j];
      }
      pm2m(ggx,ggy,ggz,ggxd,ggyd,ggzd,nmp,nfjc);
      for( j=0; j<nmp; j++ ) {
        gx[ib][j] += ggx[j];
        gy[ib][j] += ggy[j];
        gz[ib][j] += ggz[j];
      }
// M2M @ sublev 2 to npb
    } else if ( 4 <= nc3 && nc3 <= 5 ) {
      ib = ip-1;
      if( nc[2] == nc3-5 ) {
        jb = jp-1;
      } else {
        jb = jp-2;
      }
      nfjc = 7-i;
      for( j=0; j<nmp; j++ ) {
        ggxd[j] = gx[jb][j];
        ggyd[j] = gy[jb][j];
        ggzd[j] = gz[jb][j];
      }
      pm2m(ggx,ggy,ggz,ggxd,ggyd,ggzd,nmp,nfjc);
      for( j=0; j<nmp; j++ ) {
        gx[ib][j] += ggx[j];
        gy[ib][j] += ggy[j];
        gz[ib][j] += ggz[j];
      }
// M2M @ sublev 2 to npb (channel)
    } else if ( 6 <= nc3 && nc3 <= 7 ) {
      if( nc[2] == nc3-6 ) {
        ib = ip-1;
        jb = jp-1;
        nfjc = 7-i;
        for( j=0; j<nmp; j++ ) {
          ggxd[j] = gx[jb][j];
          ggyd[j] = gy[jb][j];
          ggzd[j] = gz[jb][j];
        }
        pm2m(ggx,ggy,ggz,ggxd,ggyd,ggzd,nmp,nfjc);
        for( j=0; j<nmp; j++ ) {
          gx[ib][j] += ggx[j];
          gy[ib][j] += ggy[j];
          gz[ib][j] += ggz[j];
        }
      }
// save bx for M2L @ lev 1
    } else if ( 8 <= nc3 && nc3 <= 9 ) {
      if( nc[2] == nc3-8 ) {
        ib = i+3*nbnet;
        jb = i+jp;
        for( j=0; j<nmp; j++ ) {
          gx[ib][j] += gx[jb][j];
          gy[ib][j] += gy[jb][j];
          gz[ib][j] += gz[jb][j];
        }
      }
    }
  }
}
