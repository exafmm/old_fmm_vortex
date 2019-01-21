#include "mpi.h"
#include "../misc/constants.h"

void mpibcast(int n0, int n1, float* var, int nmax) {
  int ncall,icall,iwork1,iwork2,ista,iend,nbuf,i;
  double send[mpim];

  ncall = (n1-n0+1)/mpim+1;
  for( icall=0; icall<ncall; icall++ ) {
    iwork1 = (n1-n0+1)/ncall;
    iwork2 = (n1-n0+1)%ncall;
    ista = icall*iwork1+n0+std::min(icall,iwork2);
    iend = ista+iwork1-1;
    if( iwork2 > icall ) iend++;
    nbuf = iend-ista+1;
    for( i=ista; i<=iend; i++ ) send[i-ista] = var[i];
    MPI_Bcast(send,nbuf,MPI_FLOAT,0,MPI_COMM_WORLD);
    for( i=ista; i<=iend; i++ ) var[i] = send[i-ista];
  }

  
}
