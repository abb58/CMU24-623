/*
 * brief - driver methods for support
 */
 

//------ Function Declarations ------
int get_natoms(const char* filename);
void read_xyz(const int Natoms, const char* filename);
void write_xyz(std::ofstream simFile, int config);
void dump_stats(std::ofstream enerFile, int config);

void init();
void calc_pairenergy();
void calc_virial_force(int timeunit);
void calc_kenergy();
void calc_momentum();
void calc_inst_temp_pr();
void vv_scheme();

double elapsed_time=0.0;
double dt=0.002;
double t;
int    Natoms;
double *m=NULL;
double Rcut=2.5;
double R2cut=0.0, RIcut=0.0, RI2cut=0.0, RI6cut=0.0, RI12cut=0.0;
double **r=NULL, **v=NULL, **f=NULL;
double **r_old=NULL, **v_old=NULL, **f_old=NULL;
double lx=6.8, ly=6.8, lz=6.8;
double Hlx=0.5*lx, Hly=0.5*ly, Hlz=0.5*lz;
double Vir=0.0, V=0.0;
double KB=1.0;
double T=0.0, P=0.0;
  

// Energy and force calculation variable
bool continuous=0;
double URcut=0.0, U=0.0, KE=0.0, TE=0.0;
double dx=0.0, dy=0.0, dz=0.0;
double R=0.0, r2=0.0, r6=0.0, Ir6=0.0;
double px=0.0, py=0.0, pz=0.0;
double fijx=0.0, fijy=0.0, fijz=0.0;
double F=0.0, fcut=0.0;
