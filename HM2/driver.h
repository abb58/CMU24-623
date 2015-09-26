/*
 * brief - driver methods for support
 */
 

//------ Function Declarations ------
int get_natoms(const char* filename);
void write_xyz();
void read_xyz(int natoms, const char* filename);
void dump_stats();

void init();
double calc_energy_force();
void integrate();


double dt=0.002;
double t;
double m=1;
int Natoms;
double **r=NULL, **v=NULL, **f=NULL;
double **r_old=NULL, **v_old=NULL, **f_old=NULL;

// Energy and force calculation variables
double u=0.0, ke=0.0;
double U=0.0, KE=0.0, TE=0.0;
double dx=0.0, dy=0.0, dz=0.0;
double r2=0.0, r6=0.0, Ir6=0.0;
double F=0.0;
