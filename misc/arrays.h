#include <complex>

// VM coordinates, velocity, and vortex strength

float *xi,*yi,*zi;
float *gxi,*gyi,*gzi,*vi;
float *xj,*yj,*zj;
float *gxj,*gyj,*gzj,*vj,*sj;
float *xo,*yo,*zo;
float *uo,*vo,*wo;
float *xw,*yw,*zw;
float *gxw,*gyw,*gzw;
float *uw,*vw,*ww;
float *xe,*ye,*ze;
float *gxe,*gye,*gze,*ve,*se;

// FDM velocity, pressure and index shifter

float ***u,***v,***w;
float ***um,***vm,***wm;
float ***ua,***va,***wa;
float ***uc,***vc,***wc;
float ***p,***po,***q;
float ***dpx,***dpy,***dpz;
float ***dphx,***dphy,***dphz;
double ***phi;
std::complex<double> ***phit;
int *m2,*m1,*l1,*l2;

// FDM R-K coefficients, and boundary values

double *cga,*cze,*cal;
double *vxl,*vyl,*vzl;
double *vxu,*vyu,*vzu;

// Channel FDM variables

double *x1,*x2,*x12,*dx1,*dx12;
double *y11,*y2,*y12,*dy1,*dy12;
double *z1,*z2,*z12,*dz1,*dz12;
double *cx,*cz,*cxz;
double *sta;

// Spectral method velocity

std::complex<double> ***uk,***vk,***wk;
std::complex<double> ***uko,***vko,***wko;
std::complex<double> ***uk1,***vk1,***wk1;
std::complex<double> ***uk2,***vk2,***wk2;
std::complex<double> ***uk3,***vk3,***wk3;
std::complex<double> ***uk4,***vk4,***wk4;
std::complex<double> ***fku,***fkv,***fkw;

// Turbulence statistics variables

std::complex<double> ***ux,***uy,***uz;
std::complex<double> ***vx,***vy,***vz;
std::complex<double> ***wx,***wy,***wz;

// FFTW plan and arrays

long *plan;
std::complex<float> *s1,*s2,*s3;
std::complex<float> *r1,*r2,*r3;
std::complex<double> ***uf,***vf,***wf;
std::complex<double> *u1,*v1,*w1;
std::complex<double> **u2,**v2,**w2;
std::complex<double> ***u3,***v3,***w3;
std::complex<double> ***u4,***v4,***w4;

// BEM quadratures and node values

double **ptr;
int **ntr,**nne;
double **vna,**vta,**vsa;
double **p1,**p2,**p3;
double **p4,**p5,**p6;
double **hsg,*are,*crv;
double *si,*gt,*gs;
double *un,*ut,*us;
double *xiq,*etq,*wq;
double (*dimh)[nwmax],*ust,*gts;

// FMM indexing

int *nbi,**ndi,*nei,*nfi,*nlbi;
int *nbj,**ndj,*nej,*nfj,*nlbj;
int **ndl,*nel,*nfl,*nlbl;
int **ndm,*nem,*nfm,*nlbm;
int **ndo,*neo,*nfo,**nxs,*nfn;
int *nij,*njb,*nnp,**neij,**nsij;
int **npx,**npxd,*ncnt,**ncnt2,*nfrom,*nc,*nd;

// FMM moments and translation components

std::complex<double> (*ax)[mpsym],(*axo)[mpsym],(*bx)[mpsym];
std::complex<double> (*ay)[mpsym],(*ayo)[mpsym],(*by)[mpsym];
std::complex<double> (*az)[mpsym],(*azo)[mpsym],(*bz)[mpsym];
std::complex<double> *bnm,*bth,*ynm,***dnm;
float *fac,*sr,*anm,*ank;

// PPM moments and translation components

double (*px)[nspm],(*pxo)[nspm],(*gx)[nspm];
double (*py)[nspm],(*pyo)[nspm],(*gy)[nspm];
double (*pz)[nspm],(*pzo)[nspm],(*gz)[nspm];
double *xsp,*ysp,*zsp;
double **psm,**p2g;

// FMM,PPM auxiliary variables

float *xd,*gxd,*gyd,*gzd,*vd,*sortd;
int *nsortd;
double *tfmm;

// ParMETIS variables

int *nek,*ixadj,*nxadj,*nadjncy;
int *nvwgt,*nadjwgt,*nvtxdist,*npart,*npartd;
int *irank,*mxadj,*madjncy,*mvwgt,*madjwgt;

// MDGRAPE,GPU indexing

int *ibase,*isize,*jbase,*jsize;
int *istamd,*iendmd,*jstamd,*jendmd;
int *istagp,*iendgp,*jstagp,*jendgp;
int *nicall,*njcall;
int *nvecd,*njj;

// MDGRAPE coordinates, velocity, and vortex strength

double *bxmd,*bymd,*bzmd;
double (*xmd)[3],(*ymd)[3],(*zmd)[3];
double (*pos)[3],(*pod)[3];

// GPU coordinates, velocity, and vortex strength

float *xig,*yig,*zig;
float *gxig,*gyig,*gzig,*vig;
float *xjg,*yjg,*zjg;
float *gxjg,*gyjg,*gzjg,*vjg,*sjg;
float *gxdg,*gydg,*gzdg,*vdg;

// GPU multipole & local expansions

float *arex,*arey,*arez;
float *aimx,*aimy,*aimz;
float *brex,*brey,*brez;
float *bimx,*bimy,*bimz;
float *ynmre,*ynmim;
float *dnmre,*dnmim;

// MPI indexing & variables

int *isdsp,*iscnt,*irdsp,*ircnt,*isort;
int *jsdsp,*jscnt,*jrdsp,*jrcnt,*jsort;
int *ksdsp,*kscnt,*krdsp,*krcnt;
int *lsdsp,*lscnt,*lrdsp,*lrcnt;
int *nsdsp,*nscnt,*nrdsp,*nrcnt;
int *nsend,*nrecv,*nbuf,*nbufd;
float *fsend,*frecv,*fbuf,*fbufd;
double *dsend,*drecv,*dbuf,*dbufd;
std::complex<float> *csend,*crecv,*cbuf,*cbufd;

// Matrix solver & RBF solver

double **dima;
float **work,*rbi,*rbj;
double **aa,**ac,**ae;

// Parallel quicksort variables

int *na,*nb,*natemp,*nbtemp;
int *nsta,*nend,*ibuf;
int *msta,*mend,*med,*msize;
int *lsta,*lend,*lshare;

// Spline interpolation x, y, Y, and u

double *xa,*ya,*yb,*ub;
