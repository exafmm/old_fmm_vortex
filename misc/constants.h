#define ALLOCATE_IN_MAIN
#include <cmath>
#include <complex>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sys/time.h>

const int mprint  = 0;        // 0 : don't, 1 : print memory usage
const int nwmax   = 2500;     // max of wall elements
const int npmax   = 4100000;  // max of particles
const int mpmax   = 100;      // max of FMM moments
const int mpsym   = 55;       // max of FMM moments (symmetry)
const int mpcmp   = 10;       // max of FMM moments (compressed)
const int nspm    = 216;      // max of pseudo particles
const int nebm    = 189;      // max of interacting boxes
const int nrbm    = 512;      // max of relative box positioning
const int mpim    = 500000;   // max of MPI buffer
const int mimax   = 1000000;  // max of MDGRAPE i particle memory
const int mjmax   = 50000;    // max of MDGRAPE j particle memory
const int mbmax   = 255;      // max of cells in MDGRAPE cell index
const int ngpu    = 4;        // GPUs per node
const int nimax   = 200000;   // max of GPU i particle memory
const int njmax   = 100000;   // max of GPU j particle memory
const int nblok0  = 128;      // size of GPU thread block P2P
const int nblok1  = 64;       // size of GPU thread block B2B
const int ins     = 50;       // insertion threshold for quicksort
const int itmax   = 1000;     // max of matrix interations
const int ngl     = 3;        // Gauss triangle 1,3,4,6,7,9,12,13
const float tol   = 1.0e-3;   // tolerance of residual
const float eps   = 1.0e-6;   // single precision epsilon
const float xmins = 1.0e-5;   // min of singular table
const float xmaxs = 1.0e3;    // max of singular table
const float xminn = 1.0e-3;   // min of non-singular table
const float xmaxn = 1.0e5;    // max of non-singular table
const float pi    = M_PI;

#ifdef MAIN
int ncheck;                  // check point switch
int nbmax;                   // total of FMM boxes @ lmax
int nbne;                    // non-empty FMM boxes @ lmax
int nbnp;                    // non-empty periodic boxes @ lmax
int nbnet;                   // total of nbne for all levels
int nbnes;                   // total of nbne with shear
int lmax;                    // number of FMM box divisions
int m2lrp2p;                 // ratio of kernel weight between m2l and p2p
int nsub;                    // number of boxes in subtree
int nwn;                     // number of boundary nodes
int nwe;                     // number of bounary elements
int nprocs;                  // number of processors in MPI
int myrank;                  // rank of node
int ierr;                    // error flag for MPI
int ndir;                    // environment variable Z0_DIR
int nmpi;                    // environment variable Z1_MPI
int ndev;                    // environment variable Z2_DEV
int nfmm;                    // environment variable Z3_FMM
int ngeo;                    // environment variable Z4_GEO
int nsmt;                    // environment variable Z5_SMT
int nequ;                    // environment variable Z6_EQU
int neq;                     // nequ for fmm internal use
int mem;                     // memory usage
float umem;                  // total memory usage
float stx;                   // shear rate parameter
float xmin;                  // xmin of FMM box
float ymin;                  // ymin of FMM box
float zmin;                  // zmin of FMM box
float rd;                    // length of FMM box
float dxyz;                  // reference length of vortex volume
float dx;                    // x length of vortex volume
float dy;                    // y length of vortex volume
float dz;                    // z length of vortex volume

double get_time(void) {
  struct timeval tv;
  struct timezone tz;
  gettimeofday(&tv, &tz);
  return ((double)(tv.tv_sec+tv.tv_usec*1.0e-6));
}

#else
extern int ncheck;
extern int nbmax;
extern int nbne;
extern int nbnp;
extern int nbnet;
extern int nbnes;
extern int lmax;
extern int m2lrp2p;
extern int nsub;
extern int nwn;
extern int nwe;
extern int nprocs;
extern int myrank;
extern int ierr;
extern int ndir;
extern int nmpi;
extern int ndev;
extern int nfmm;
extern int ngeo;
extern int nsmt;
extern int nequ;
extern int neq;
extern int mem;
extern float umem;
extern float stx;
extern float xmin;
extern float ymin;
extern float zmin;
extern float rd;
extern float dxyz;
extern float dx;
extern float dy;
extern float dz;

extern double get_time(void);
#endif

template<typename write_type>
void binary_write(std::fstream& fid, write_type *data, const int n) {
  int ibuf,i,j;
  const int nbuf=200000/sizeof(write_type);
  write_type buffer[nbuf];

  for( i=0; i<(n+nbuf-1)/nbuf; i++ ) {
    ibuf = std::min(nbuf,n-i*nbuf);
    for( j=0; j<ibuf; j++ ) buffer[j] = data[j+i*nbuf];
    fid.write((char *)(&buffer),sizeof(write_type)*ibuf);
  }
}

template<typename read_type>
void binary_read(std::fstream& fid, read_type *data, const int n) {
  int ibuf,i,j;
  const int nbuf=200000/sizeof(read_type);
  read_type buffer[nbuf];

  for( i=0; i<(n+nbuf-1)/nbuf; i++ ) {
    ibuf = std::min(nbuf,n-i*nbuf);
    fid.read((char *)(&buffer),sizeof(read_type)*ibuf);
    for( j=0; j<ibuf; j++ ) data[j+i*nbuf] = buffer[j];
  }
}
