/*
 * brief - driver methods for support
 */
 

//------ Function Declarations ------
int get_natoms(const char* filename);
void write_xyz();
void read_xyz(int natoms, const char* filename);


void init();
void calc_energy_force();
void calc_kenergy();
void vv_scheme();

double elapsed_time=0.0;
double dt=0.002;
double t;
double m=1;
int Natoms;
double **r=NULL, **v=NULL, **f=NULL;
double **r_old=NULL, **v_old=NULL, **f_old=NULL;

// Energy and force calculation variables
double U=0.0, KE=0.0, TE=0.0;
double dx=0.0, dy=0.0, dz=0.0;
double r2=0.0, r6=0.0, Ir6=0.0;
double px=0.0, py=0.0, pz=0.0;
double F=0.0;
