#include "../misc/constants.h"

extern float *xi,*yi,*zi,*gxi,*gyi,*gzi;
extern int **ndi,*nfj,*nc,**npx,**neij,*nij,*njb;
extern std::complex<double> (*bx)[mpsym],(*by)[mpsym],(*bz)[mpsym];
extern float *fac,*gxd,*gyd,*gzd;
extern int *jbase,*jsize,*istagp,*iendgp,**jstagp,**jendgp,*njcall,*nvecd,*njj;
extern float *xig,*yig,*zig,*gxig,*gyig,*gzig,*vig;
extern float *brex,*brey,*brez,*bimx,*bimy,*bimz;

extern void boxc(int, int, int*);
extern "C" void m2pgpu_(int*, double*, double*, double*, double*, double*, double*,
        float*, float*, float*, float*, float*, float*, float*,
        float*, float*, float*, float*, float*, float*, float*);

void stm2p(int nmp, int mp, int lbi, int lbj, int lev, int ipb, double rb) {
  int idev,mblok,ni,nj,nicall,jc,jj,ii,njd,ij,icall,jcall,iblok,jjd,jb,j,ibase,isize,is,i,ijc,jjdd;
  double op,rbd,epsd,xmind,ymind,zmind;

  idev = myrank%ngpu;
  mblok = 4*nebm+1;

  ni = 0;
  nj = 0;
  nicall = 0;
  istagp[0] = 0;
  jc = 0;
  for( jj=0; jj<nbnp; jj++ ) njj[jj] = 0;
  for( ii=0; ii<lbi; ii++ ) {
    if( nij[ii] != 0 ) {
      njd = 0;
      jc = 0;
      jstagp[0][ii] = 0;
      for( ij=0; ij<nij[ii]; ij++ ) {
        jj = neij[ij][ii];
        if( njj[jj] == 0 ) {
          nj += nmp;
          njj[jj] = 1;
        }
        njd += nmp;
        if( njd > njmax ) {
          jendgp[jc][ii] = ij-1;
          jc++;
          jstagp[jc][ii] = ij;
          njd = nmp;
        }
      }
      jendgp[jc][ii] = nij[ii]-1;
      ni += ((ndi[1][ii]-ndi[0][ii]+nblok1)/nblok1+1)*nblok1;
      if( jc != 0 ) {
        if( ii > istagp[nicall] ) {
          njcall[nicall] = 1;
          iendgp[nicall] = ii-1;
          nicall++;
          istagp[nicall] = ii;
          for( jj=0; jj<nbnp; jj++ ) njj[jj] = 0;
        }
        if( ii != lbi ) {
          njcall[nicall] = jc+1;
          iendgp[nicall] = ii;
          nicall++;
          istagp[nicall] = ii+1;
          for( jj=0; jj<nbnp; jj++ ) njj[jj] = 0;
          ni = 0;
          nj = 0;
        }
      } else if ( ni > nimax || nj > njmax ) {
        njcall[nicall] = jc+1;
        iendgp[nicall] = ii-1;
        nicall++;
        istagp[nicall] = ii;
        for( jj=0; jj<nbnp; jj++ ) njj[jj] = 0;
        ni = ((ndi[1][ii]-ndi[0][ii]+nblok1)/nblok1+1)*nblok1;
        nj = 0;
        for( ij=0; ij<nij[ii]; ij++ ) {
          jj = neij[ij][ii];
          nj += nmp;
          njj[jj] = 1;
        }
      }
    }
  }
  njcall[nicall] = jc+1;
  iendgp[nicall] = lbi-1;
  nicall++;

  for( icall=0; icall<nicall; icall++ ) {
    for( jcall=0; jcall<njcall[icall]; jcall++ ) {
      iblok = 0;
      jc = 0;
      jjd = 0;
      op = 0;
      for( jj=0; jj<nbnp; jj++ ) njj[jj] = 0;
      for( ii=istagp[icall]; ii<=iendgp[icall]; ii++ ) {
        if( nij[ii] != 0 ) {
          for( ij=jstagp[jcall][ii]; ij<=jendgp[jcall][ii]; ij++ ) {
            jj = neij[ij][ii];
            jb = njb[jj];
            if( njj[jj] == 0 ) {
              jbase[jjd] = jc;
              for( j=0; j<nmp; j++ ) {
                brex[jc] = std::real(bx[jb][j]);
                brey[jc] = std::real(by[jb][j]);
                brez[jc] = std::real(bz[jb][j]);
                bimx[jc] = std::imag(bx[jb][j]);
                bimy[jc] = std::imag(by[jb][j]);
                bimz[jc] = std::imag(bz[jb][j]);
                jc++;
              }
              jsize[jjd] = jc-jbase[jjd];
              jjd++;
              njj[jj] = jjd;
            }
          }
          ibase = ndi[0][ii];
          isize = ndi[1][ii]-ibase+1;
          for( is=0; is<isize; is+=nblok1 ) {
            for( i=0; i<std::min(isize-is,nblok1); i++ ) {
              xig[iblok*nblok1+i] = xi[ibase+is+i];
              yig[iblok*nblok1+i] = yi[ibase+is+i];
              zig[iblok*nblok1+i] = zi[ibase+is+i];
              gxig[iblok*nblok1+i] = gxi[ibase+is+i];
              gyig[iblok*nblok1+i] = gyi[ibase+is+i];
              gzig[iblok*nblok1+i] = gzi[ibase+is+i];
            }
            for( i=isize-is; i<nblok1; i++ ) {
              xig[iblok*nblok1+i] = 0;
              yig[iblok*nblok1+i] = 0;
              zig[iblok*nblok1+i] = 0;
              gxig[iblok*nblok1+i] = 0;
              gyig[iblok*nblok1+i] = 0;
              gzig[iblok*nblok1+i] = 0;
            }
            nvecd[iblok*mblok+10] = jendgp[jcall][ii]-jstagp[jcall][ii]+1;
            ijc = 0;
            for( ij=jstagp[jcall][ii]; ij<=jendgp[jcall][ii]; ij++ ) {
              jj = neij[ij][ii];
              if( njj[jj] != 0 ) {
                jjdd = njj[jj]-1;
                ijc++;
                boxc(nfj[jj],3,nc);
                nvecd[iblok*mblok+4*ijc+7] = jbase[jjdd];
                nvecd[iblok*mblok+4*ijc+8] = nc[0]+npx[0][jj]*pow(2,lev);
                nvecd[iblok*mblok+4*ijc+9] = nc[1]+npx[1][jj]*pow(2,lev);
                nvecd[iblok*mblok+4*ijc+10] = nc[2]+npx[2][jj]*pow(2,lev);
                op += (double) nblok1*jsize[jjdd];
              }
            }
            iblok++;
          }
        }
      }
      nvecd[0] = idev;
      nvecd[1] = iblok;
      nvecd[2] = mblok;
      nvecd[3] = jc;
      nvecd[4] = 2;
      nvecd[5] = myrank;
      nvecd[6] = mp;
      nvecd[7] = ipb;
      rbd = rb;
      epsd = eps;
      xmind = xmin;
      ymind = ymin;
      zmind = zmin;
      m2pgpu_(nvecd,&op,&rbd,&epsd,&xmind,&ymind,&zmind,xig,yig,zig,gxig,gyig,gzig,vig,
              brex,brey,brez,bimx,bimy,bimz,fac);
      iblok = 0;
      for( ii=istagp[icall]; ii<=iendgp[icall]; ii++ ) {
        if( nij[ii] != 0 ) {
          ibase = ndi[0][ii];
          isize = ndi[1][ii]-ibase+1;
          for( is=0; is<isize; is+=nblok1 ) {
            for( i=0; i<std::min(isize-is,nblok1); i++ ) {
              gxd[ibase+is+i] += gxig[iblok*nblok1+i];
              gyd[ibase+is+i] += gyig[iblok*nblok1+i];
              gzd[ibase+is+i] += gzig[iblok*nblok1+i];
            }
            iblok++;
          }
        }
      }
    }
  }
}
