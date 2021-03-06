#ifdef MAIN
int    ndv = 2;      // number of wall element divisions
int    nge = 5;      // geometry number
int    nin = 1;      // 1 stop 2 erase
int    nre = 1;      // 1 don't 2 remesh
int    nst = 2;      // 1 don't 2 calculate statistics
int    mp  = 10;     // order of multipole moments
int    npb = 0;      // exponent of periodic images
double rw  = 1.0;    // radius of sphere
double rwb = 0.75;   // radius of boundary
double dt  = 0.1;    // time step size
double hs  = 0.05;   // height of initial element
double uin = 1.0;    // inlet velocity
double vin = 0.0;    // inlet velocity
double win = 0.0;    // inlet velocity
double vis = 0.001;  // kinematic viscosity
double cgp = 10;     // strength of intial element
#else
extern int    ndv = 2;      // number of wall element divisions
extern int    nge = 5;      // geometry number
extern int    nin = 1;      // 1 stop 2 erase
extern int    nre = 1;      // 1 don't 2 remesh
extern int    nst = 2;      // 1 don't 2 calculate statistics
extern int    mp  = 10;     // order of multipole moments
extern int    npb = 0;      // exponent of periodic images
extern double rw  = 1.0;    // radius of sphere
extern double rwb = 0.75;   // radius of boundary
extern double dt  = 0.1;    // time step size
extern double hs  = 0.05;   // height of initial element
extern double uin = 1.0;    // inlet velocity
extern double vin = 0.0;    // inlet velocity
extern double win = 0.0;    // inlet velocity
extern double vis = 0.001;  // kinematic viscosity
extern double cgp = 10;     // strength of intial element
#endif
