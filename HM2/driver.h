/*
 * brief - driver methods for support
 */
 

//------ Function Declarations ------
int get_natoms(const char* filename);
void write_xyz();
void read_xyz(int natoms, const char* filename);

void print_energy();
void init();
double tot_energy();
double pot_energy();
double kin_energy();
void integrate();

double dt=0.002;
double t;
double m=1;
int Natoms;
double **r=NULL, **v=NULL, **f=NULL;
double **r_old=NULL, **v_old=NULL, **f_old=NULL;

