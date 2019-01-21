const int nge   = 0;       // geometry number
const int mp    = 10;      // order of multipole moments
const double dt  = 1e-2;   // dummy time step size
const double vis = 2e-2;   // kinematic viscosity
#ifdef MAIN
int npb,ngg;
#else
extern int npb,ngg;
#endif
