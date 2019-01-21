#include "../misc/constants.h"

extern int *nc;
extern double (*px)[nspm],(*py)[nspm],(*pz)[nspm],*xsp,*ysp,*zsp,**psm,**p2g;

extern void boxc(int, int, int*);
extern void inv(double**, int);

void trans(int& nmp, int mp) {
  int i,j,k,n,ix,iy,iz,ni,njk,nbyte,msp[10];
  double xd,yd,zd,rspd,xjjc,yjjc,zjjc,rh,cosg,psmd,dxij,dyij,dzij,rij,pn[mpmax];
  std::fstream fid;
  fid.open("../../dat/sphericalt",std::ios::in);

  msp[0] = 4;
  msp[1] = 14;
  msp[2] = 26;
  msp[3] = 36;
  msp[4] = 60;
  msp[5] = 84;
  msp[6] = 108;
  msp[7] = 144;
  msp[8] = 180;
  msp[9] = 216;

  nmp = 0;
  for( i=0; i<mp-1; i++ ) nmp += msp[i];
  for( i=0; i<nmp; i++ ) {
    fid >> xd;
    fid >> yd;
    fid >> zd;
  }
  for( i=0; i<msp[mp-1]; i++ ) {
    fid >> xsp[i];
    fid >> ysp[i];
    fid >> zsp[i];
  }
  nmp = msp[mp-1];

  fid.close();

  rspd = sqrt(3.0)*0.5;
  pn[0] = 0;
  pn[1] = 1;
  for( iz=0; iz<2; iz++ ) {
    for( ix=0; ix<2; ix++ ) {
      for( iy=0; iy<2; iy++ ) {
        ni = iz*4+ix*2+iy;
        for( j=0; j<nmp; j++ ) {
          xjjc = rspd*xsp[j]+0.5*pow(-1,ix+1);
          yjjc = rspd*ysp[j]+0.5*pow(-1,iy+1);
          zjjc = rspd*zsp[j]+0.5*pow(-1,iz+1);
          rh = sqrt(xjjc*xjjc+yjjc*yjjc+zjjc*zjjc)+eps;
          for( k=0; k<nmp; k++ ) {
            njk = j*nmp+k;
            cosg = (xjjc*xsp[k]+yjjc*ysp[k]+zjjc*zsp[k])/rh;
            psmd = 1.0/nmp;
            for( n=1; n<mp; n++ ) {
              pn[n+1] = (2*n-1.0)/n*cosg*pn[n]-(n-1.0)/n*pn[n-1];
              psmd += (2*n+1.0)/nmp*pow(rh*0.5/rspd,n)*pn[n+1];
            }
            psm[ni][njk] = psmd;
          }
        }
      }
    }
  }

  for( i=0; i<nmp; i++ ) {
    for( j=0; j<nmp; j++ ) {
      dxij = 2*xsp[i]-xsp[j];
      dyij = 2*ysp[i]-ysp[j];
      dzij = 2*zsp[i]-zsp[j];
      rij = sqrt(dxij*dxij+dyij*dyij+dzij*dzij)+eps;
      p2g[i][j] = 0.25/pi/rij;
    }
  }
  inv(p2g,nmp);

  for( j=0; j<nbnes; j++ ) {
    for( i=0; i<nmp; i++ ) {
      px[j][i] = 0;
      py[j][i] = 0;
      pz[j][i] = 0;
    }
  }
}
