#include "../misc/parameters.h"
#include "../misc/constants.h"

extern float *xi,*yi,*zi,*gxi,*gyi,*gzi,*vi;
extern float *xj,*yj,*zj,*gxj,*gyj,*gzj,*vj,*sj;
extern int **ndi,*nei,**ndj,*nej,*nij,**neij;
extern float *gxd,*gyd,*gzd,*vd;
extern int *jbase,*jsize,*istagp,*iendgp,**jstagp,**jendgp,*njcall,*nvecd,*njj;
extern float *xig,*yig,*zig,*gxig,*gyig,*gzig,*vig;
extern float *xjg,*yjg,*zjg,*gxjg,*gyjg,*gzjg,*vjg,*sjg;
 extern float *gxdg,*gydg,*gzdg,*vdg;

extern "C" void p2pgpu_(int*, double*, double*, double*, double*, double*, double*, double*,
        float*, float*, float*, float*, float*, float*, float*,
        float*, float*, float*, float*, float*, float*, float*, float*,
        float*, float*, float*, float*);

void p2pgp(int lbi, int neqd) {
  int num,idev,mblok,ni,nj,nicall,jc,jj,ii,njd,ij,icall,jcall,iblok,jjd,j,ibase,isize,is,i,ijc,jjdd;
  double op,visd,epsd,dxd,dyd,dzd,dxyzd;

  num = nprocs/pow(8,lmax);
  if( num == 0 ) num = 1;
  idev = myrank%ngpu;
  mblok = 60*num;

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
          nj += ndj[1][jj]-ndj[0][jj]+1;
          njj[jj] = 1;
        }
        njd += ndj[1][jj]-ndj[0][jj]+1;
        if( njd > njmax ) {
          jendgp[jc][ii] = ij-1;
          jc++;
          jstagp[jc][ii] = ij;
          njd = ndj[1][jj]-ndj[0][jj]+1;
        }
      }
      jendgp[jc][ii] = nij[ii]-1;
      ni += std::max(ndi[1][ii]-ndi[0][ii]+1,nblok0);
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
        ni = std::max(ndi[1][ii]-ndi[0][ii]+1,nblok0);
        nj = 0;
        for( ij=0; ij<nij[ii]; ij++ ) {
          jj = neij[ij][ii];
          nj += ndj[1][jj]-ndj[0][jj]+1;
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
            if( njj[jj] == 0 ) {
              jbase[jjd] = jc;
              for( j=ndj[0][jj]; j<=ndj[1][jj]; j++ ) {
                xjg[jc] = xj[j];
                yjg[jc] = yj[j];
                zjg[jc] = zj[j];
                gxjg[jc] = gxj[j];
                gyjg[jc] = gyj[j];
                gzjg[jc] = gzj[j];
                vjg[jc] = vj[j];
                sjg[jc] = sj[j];
                jc++;
              }
              jsize[jjd] = jc-jbase[jjd];
              jjd++;
              njj[jj] = jjd;
            }
          }
          ibase = ndi[0][ii];
          isize = ndi[1][ii]-ibase+1;
          for( is=0; is<isize; is+=nblok0 ) {
            for( i=0; i<std::min(isize-is,nblok0); i++ ) {
              xig[iblok*nblok0+i] = xi[ibase+is+i];
              yig[iblok*nblok0+i] = yi[ibase+is+i];
              zig[iblok*nblok0+i] = zi[ibase+is+i];
              gxig[iblok*nblok0+i] = gxi[ibase+is+i];
              gyig[iblok*nblok0+i] = gyi[ibase+is+i];
              gzig[iblok*nblok0+i] = gzi[ibase+is+i];
              vig[iblok*nblok0+i] = vi[ibase+is+i];
            }
            for( i=isize-is; i<nblok0; i++ ) {
              xig[iblok*nblok0+i] = 0;
              yig[iblok*nblok0+i] = 0;
              zig[iblok*nblok0+i] = 0;
              gxig[iblok*nblok0+i] = 0;
              gyig[iblok*nblok0+i] = 0;
              gzig[iblok*nblok0+i] = 0;
              vig[iblok*nblok0+i] = 0;
            }
            nvecd[iblok*mblok+10] = jendgp[jcall][ii]-jstagp[jcall][ii]+1;
            ijc = 0;
            for( ij=jstagp[jcall][ii]; ij<=jendgp[jcall][ii]; ij++ ) {
              jj = neij[ij][ii];
              if( njj[jj] != 0 ) {
                jjdd = njj[jj]-1;
                ijc++;
                nvecd[iblok*mblok+2*ijc+9] = jbase[jjdd];
                nvecd[iblok*mblok+2*ijc+10] = jsize[jjdd];
                op += (double) nblok0*jsize[jjdd];
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
      nvecd[4] = neqd;
      nvecd[5] = myrank;
      visd = vis;
      epsd = eps;
      dxd = dx;
      dyd = dy;
      dzd = dz;
      dxyzd = dxyz;
      p2pgpu_(nvecd,&op,&visd,&epsd,&dxd,&dyd,&dzd,&dxyzd,xig,yig,zig,gxig,gyig,gzig,vig,
              xjg,yjg,zjg,gxjg,gyjg,gzjg,vjg,sjg,gxdg,gydg,gzdg,vdg);
      iblok = 0;
      for( ii=istagp[icall]; ii<=iendgp[icall]; ii++ ) {
        if( nij[ii] != 0 ) {
          ibase = ndi[0][ii];
          isize = ndi[1][ii]-ibase+1;
          for( is=0; is<isize; is+=nblok0 ) {
            for( i=0; i<std::min(isize-is,nblok0); i++ ) {
              gxd[ibase+is+i] += gxdg[iblok*nblok0+i];
              gyd[ibase+is+i] += gydg[iblok*nblok0+i];
              gzd[ibase+is+i] += gzdg[iblok*nblok0+i];
              vd[ibase+is+i] += vdg[iblok*nblok0+i];
            }
            iblok++;
          }
        }
      }
    }
  }
}
