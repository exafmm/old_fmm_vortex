#include "../misc/constants.h"

extern int *nfi,**ndj,*nfj,**neij,*nij,*nek,*ixadj,*nadjncy,*nadjwgt,*npart;

void setedge(int lbi, int neib) {
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
      if( iv != jv && npart[jv] == 0 ) {
        nadjncy[ixadj[iv]] = jv;
        nadjncy[ixadj[jv]] = iv;
        if( neib == 2 ) {
          nadjwgt[ixadj[iv]] += (ndj[1][jj]-ndj[0][jj]+1);
          nadjwgt[ixadj[jv]] += (ndj[1][jj]-ndj[0][jj]+1);
        } else if( neib == 4 ) {
          nadjwgt[ixadj[iv]] += mpsym*3/8;
          nadjwgt[ixadj[jv]] += mpsym*3/8;
        }
        ixadj[iv]++;
        ixadj[jv]++;
        npart[jv] = 1;
      }
    }
  }
}
