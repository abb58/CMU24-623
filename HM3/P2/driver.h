/*
 * brief - driver methods for support
 */
 

//------ Function Declarations ------
int get_natoms(const char* filename);
void write_xyz();
void read_xyz(int natoms, const char* filename);


void init();
void calc_pairenergy();
void calc_force();
void calc_kenergy();
void vv_scheme();

double elapsed_time=0.0;
double dt=0.002;
double t;
int    Natoms;
double *m=NULL;
double rcut=2.5;
double R2cut=0.0, RIcut=0.0, RI2cut=0.0, RI6cut=0.0, RI12cut=0.0;
double **r=NULL, **v=NULL, **f=NULL;
double **r_old=NULL, **v_old=NULL, **f_old=NULL;
double lx=6.8, ly=6.8, lz=6.8;
double Hlx=0.5*lx, Hly=0.5*ly, Hlz=0.5*lz;

// Energy and force calculation variables
double U=0.0, KE=0.0, TE=0.0;
double dx=0.0, dy=0.0, dz=0.0;
double R=0.0, r2=0.0, r6=0.0, Ir6=0.0;
double px=0.0, py=0.0, pz=0.0;
double F=0.0, fcut=0.0;
