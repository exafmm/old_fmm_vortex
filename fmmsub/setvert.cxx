#include "../misc/constants.h"

extern int **ndi,*nfi,**ndj,*nfj,**neij,*nij,*nek,*nxadj,*nvwgt,*npart;

void setvert(int lbi, int neib) {
  int ic,ii,iv,i,ij,jj,jv;

  ic = -1;
  for( ii=0; ii<lbi; ii++ ) {
    iv = (nek[nfi[ii]])/nsub;
    if( ic != iv ) {
      ic = iv;
      for( i=0; i<nbne; i++ ) npart[i] = 0;
    }
    for( ij=0; ij<nij[ii]; ij++ ) {
      jj = neij[ij][ii];
      jv = nek[nfj[jj]]/nsub;
      if( neib == 2 ) {
        nvwgt[iv] += (ndi[1][ii]-ndi[0][ii]+1)*(ndj[1][jj]-ndj[0][jj]+1);
      } else if( neib == 4 ) {
        nvwgt[iv] += m2lrp2p;
      }
      if( iv != jv && npart[jv] == 0 ) {
        nxadj[iv+1]++;
        nxadj[jv+1]++;
        npart[jv] = 1;
      }
    }
  }
}
