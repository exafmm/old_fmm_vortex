#include "mpi.h"
#include "../misc/constants.h"

extern int *nbuf,*nbufd;
extern float *fbuf,*fbufd;
extern double *dbuf,*dbufd;
extern std::complex<float> *cbuf,*cbufd;

extern void memoryuse();
extern void memoryfree();
extern void mpisendi(int*, int, int, int);
extern void mpisendrecvi(int, int, int, int);
extern void mpirecvi(int*, int, int, int);
extern void mpisendf(float*, int, int, int);
extern void mpisendrecvf(int, int, int, int);
extern void mpirecvf(float*, int, int, int);
extern void mpisendd(double*, int, int, int);
extern void mpisendrecvd(int, int, int, int);
extern void mpirecvd(double*, int, int, int);
extern void mpisendc(std::complex<float>*, int, int, int);
extern void mpisendrecvc(int, int, int, int);
extern void mpirecvc(std::complex<float>*, int, int, int);

// mpi_alltoallv for integer
void mpialltoallvi(int* nvar, int* nscnt, int* nsdsp, int* nvard, int* nrcnt, int* nrdsp, int nmax) {
  int lsend,lrecv,isize,jsize,ista,iend,jsta,jend,irank,i;
  MPI_Status istatus;
  MPI_Request ireq1,ireq2;
  
  nbuf = new int [npmax];
  mem = npmax*4;
  memoryuse();

  MPI_Barrier(MPI_COMM_WORLD);
  if( nmpi == 1 ) {
    MPI_Alltoallv(nvar,nscnt,nsdsp,MPI_INT,nvard,nrcnt,nrdsp,MPI_INT,MPI_COMM_WORLD);
  } else if( nmpi == 2 ) {
    for( lsend=0; lsend<nprocs; lsend++ ) {
      for( lrecv=0; lrecv<nprocs; lrecv++ ) {
        isize = 0;
        jsize = 0;
        if( myrank == lsend ) {
          ista = nsdsp[lrecv];
          isize = nscnt[lrecv];
          iend = ista+isize-1;
        }
        if( myrank == lrecv ) {
          jsta = nrdsp[lsend];
          isize = nrcnt[lsend];
          jend = jsta+isize-1;
        }
        if( isize > 0 ) {
          mpisendi(nvar,ista,iend,lsend);
          mpisendrecvi(isize,lsend,lrecv,1);
          mpirecvi(nvard,jsta,jend,lrecv);
        }
      }
    }
  } else if( nmpi == 3 ) {
    nbufd = new int [npmax];
    mem = npmax*4;
    memoryuse();
    for( irank=0; irank<nprocs; irank++ ) {
      lrecv = (myrank+irank)%nprocs;
      lsend = (myrank-irank+nprocs)%nprocs;
      ista = nsdsp[lrecv];
      isize = nscnt[lrecv];
      iend = ista+isize-1;
      jsta = nrdsp[lsend];
      jsize = nrcnt[lsend];
      jend = jsta+jsize-1;
      for( i=ista; i<=iend; i++ ) nbuf[i-ista] = nvar[i];
      if( irank != 0 ) {
        MPI_Isend(nbuf,isize,MPI_INT,lrecv,1,MPI_COMM_WORLD,&ireq1);
        MPI_Irecv(nbufd,jsize,MPI_INT,lsend,1,MPI_COMM_WORLD,&ireq2);
        MPI_Wait(&ireq1,&istatus);
        MPI_Wait(&ireq2,&istatus);
      } else {
        for( i=0; i<isize; i++ ) nbufd[i] = nbuf[i];
      }
      for( i=jsta; i<=jend; i++ ) nvard[i] = nbufd[i-jsta];
    }
    delete[] nbufd;
    mem = npmax*4;
    memoryfree();
  }
  delete[] nbuf;
  mem = npmax*4;
  memoryfree();
}

// mpi_alltoallv for float
void mpialltoallvf(float* fvar, int* nscnt, int* nsdsp, float* fvard, int* nrcnt, int* nrdsp, int nmax) {
  int lsend,lrecv,isize,jsize,ista,iend,jsta,jend,irank,i;
  MPI_Status istatus;
  MPI_Request ireq1,ireq2;

  fbuf = new float [npmax];
  mem = npmax*4;
  memoryuse();

  MPI_Barrier(MPI_COMM_WORLD);
  if( nmpi == 1 ) {
    MPI_Alltoallv(fvar,nscnt,nsdsp,MPI_FLOAT,fvard,nrcnt,nrdsp,MPI_FLOAT,MPI_COMM_WORLD);
  } else if( nmpi == 2 ) {
    for( lsend=0; lsend<nprocs; lsend++ ) {
      for( lrecv=0; lrecv<nprocs; lrecv++ ) {
        isize = 0;
        jsize = 0;
        if( myrank == lsend ) {
          ista = nsdsp[lrecv];
          isize = nscnt[lrecv];
          iend = ista+isize-1;
        }
        if( myrank == lrecv ) {
          jsta = nrdsp[lsend];
          isize = nrcnt[lsend];
          jend = jsta+isize-1;
        }
        if( isize > 0 ) {
          mpisendf(fvar,ista,iend,lsend);
          mpisendrecvf(isize,lsend,lrecv,1);
          mpirecvf(fvard,jsta,jend,lrecv);
        }
      }
    }
  } else if( nmpi == 3 ) {
    fbufd = new float [npmax];
    mem = npmax*4;
    memoryuse();
    for( irank=0; irank<nprocs; irank++ ) {
      lrecv = (myrank+irank)%nprocs;
      lsend = (myrank-irank+nprocs)%nprocs;
      ista = nsdsp[lrecv];
      isize = nscnt[lrecv];
      iend = ista+isize-1;
      jsta = nrdsp[lsend];
      jsize = nrcnt[lsend];
      jend = jsta+jsize-1;
      for( i=ista; i<=iend; i++ ) fbuf[i-ista] = fvar[i];
      if( irank != 0 ) {
        MPI_Isend(fbuf,isize,MPI_FLOAT,lrecv,1,MPI_COMM_WORLD,&ireq1);
        MPI_Irecv(fbufd,jsize,MPI_FLOAT,lsend,1,MPI_COMM_WORLD,&ireq2);
        MPI_Wait(&ireq1,&istatus);
        MPI_Wait(&ireq2,&istatus);
      } else {
        for( i=0; i<isize; i++ ) fbufd[i] = fbuf[i];
      }
      for( i=jsta; i<=jend; i++ ) fvard[i] = fbufd[i-jsta];
    } 
    delete[] fbufd;
    mem = npmax*4;
    memoryfree();
  }
  delete[] fbuf;
  mem = npmax*4;
  memoryfree();
}

// mpi_alltoallv for double
void mpialltoallvd(double* dvar, int* nscnt, int* nsdsp, double* dvard, int* nrcnt, int* nrdsp, int nmax) {
  int lsend,lrecv,isize,jsize,ista,iend,jsta,jend,irank,i;
  MPI_Status istatus;
  MPI_Request ireq1,ireq2;

  dbuf = new double [npmax];
  mem = npmax*8;
  memoryuse();

  MPI_Barrier(MPI_COMM_WORLD);
  if( nmpi == 1 ) {
    MPI_Alltoallv(dvar,nscnt,nsdsp,MPI_DOUBLE,dvard,nrcnt,nrdsp,MPI_DOUBLE,MPI_COMM_WORLD);
  } else if( nmpi == 2 ) {
    for( lsend=0; lsend<nprocs; lsend++ ) {
      for( lrecv=0; lrecv<nprocs; lrecv++ ) {
        isize = 0;
        jsize = 0;
        if( myrank == lsend ) {
          ista = nsdsp[lrecv];
          isize = nscnt[lrecv];
          iend = ista+isize-1;
        }
        if( myrank == lrecv ) {
          jsta = nrdsp[lsend];
          isize = nrcnt[lsend];
          jend = jsta+isize-1;
        }
        if( isize > 0 ) {
          mpisendd(dvar,ista,iend,lsend);
          mpisendrecvd(isize,lsend,lrecv,1);
          mpirecvd(dvard,jsta,jend,lrecv);
        }
      }
    }
  } else if( nmpi == 3 ) {
    dbufd = new double [npmax];
    mem = npmax*8;
    memoryuse();
    for( irank=0; irank<nprocs; irank++ ) {
      lrecv = (myrank+irank)%nprocs;
      lsend = (myrank-irank+nprocs)%nprocs;
      ista = nsdsp[lrecv];
      isize = nscnt[lrecv];
      iend = ista+isize-1;
      jsta = nrdsp[lsend];
      jsize = nrcnt[lsend];
      jend = jsta+jsize-1;
      for( i=ista; i<=iend; i++ ) dbuf[i-ista] = dvar[i];
      if( irank != 0 ) {
        MPI_Isend(dbuf,isize,MPI_DOUBLE,lrecv,1,MPI_COMM_WORLD,&ireq1);
        MPI_Irecv(dbufd,jsize,MPI_DOUBLE,lsend,1,MPI_COMM_WORLD,&ireq2);
        MPI_Wait(&ireq1,&istatus);
        MPI_Wait(&ireq2,&istatus);
      } else {
        for( i=0; i<isize; i++ ) dbufd[i] = dbuf[i];
      }
      for( i=jsta; i<=jend; i++ ) dvard[i] = dbufd[i-jsta];
    }
    delete[] dbufd;
    mem = npmax*8;
    memoryfree();
  }
  delete[] dbuf;
  mem = npmax*8;
  memoryfree();
}

// mpi_alltoallv for complex
void mpialltoallvc(std::complex<float>* cvar, int* nscnt, int* nsdsp, std::complex<float>* cvard, int* nrcnt, int* nrdsp, int nmax) {
  int lsend,lrecv,isize,jsize,ista,iend,jsta,jend,irank,i;
  MPI_Status istatus;
  MPI_Request ireq1,ireq2;

  cbuf = new std::complex<float> [npmax];
  mem = npmax*8;
  memoryuse();

  MPI_Barrier(MPI_COMM_WORLD);
  if( nmpi == 1 ) {
    MPI_Alltoallv(cvar,nscnt,nsdsp,MPI_COMPLEX,cvard,nrcnt,nrdsp,MPI_COMPLEX,MPI_COMM_WORLD);
  } else if( nmpi == 2 ) {
    for( lsend=0; lsend<nprocs; lsend++ ) {
      for( lrecv=0; lrecv<nprocs; lrecv++ ) {
        isize = 0;
        jsize = 0;
        if( myrank == lsend ) {
          ista = nsdsp[lrecv];
          isize = nscnt[lrecv];
          iend = ista+isize-1;
        }
        if( myrank == lrecv ) {
          jsta = nrdsp[lsend];
          isize = nrcnt[lsend];
          jend = jsta+isize-1;
        }
        if( isize > 0 ) {
          mpisendc(cvar,ista,iend,lsend);
          mpisendrecvc(isize,lsend,lrecv,1);
          mpirecvc(cvard,jsta,jend,lrecv);
        }
      }
    }
  } else if( nmpi == 3 ) {
    cbufd = new std::complex<float> [npmax];
    mem = npmax*8;
    memoryuse();
    for( irank=0; irank<nprocs; irank++ ) {
      lrecv = (myrank+irank)%nprocs;
      lsend = (myrank-irank+nprocs)%nprocs;
      ista = nsdsp[lrecv];
      isize = nscnt[lrecv];
      iend = ista+isize-1;
      jsta = nrdsp[lsend];
      jsize = nrcnt[lsend];
      jend = jsta+jsize-1;
      for( i=ista; i<=iend; i++ ) cbuf[i-ista] = cvar[i];
      if( irank != 0 ) {
        MPI_Isend(cbuf,isize,MPI_COMPLEX,lrecv,1,MPI_COMM_WORLD,&ireq1);
        MPI_Irecv(cbufd,jsize,MPI_COMPLEX,lsend,1,MPI_COMM_WORLD,&ireq2);
        MPI_Wait(&ireq1,&istatus);
        MPI_Wait(&ireq2,&istatus);
      } else {
        for( i=0; i<isize; i++ ) cbufd[i] = cbuf[i];
      }
      for( i=jsta; i<=jend; i++ ) cvard[i] = cbufd[i-jsta];
    }
    delete[] cbufd;
    mem = npmax*8;
    memoryfree();
  }
  delete[] cbuf;
  mem = npmax*8;
  memoryfree();
}
