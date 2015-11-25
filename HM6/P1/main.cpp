/*
 * Metropolis Monte-Carlo simulation of lennard-jones nano-
 * particles.
 *
 * Abhishek Bagusetty
 *
 * Written for CMU24-623 Molecular Simulation of Materials,
 * Fall 2015-16
 *
 * Building instructions :
 *    - make
 *    - make clean
 *
 * (c) 2015-16
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <iomanip>
#include <stdexcept>

// Function definitions
double compute_U();
void compute_delU();
void compute_trial_params();
void calc_virial_force();
void calc_inst_pressure();


double U_post=0.0, U_pre=0.0;

double elapsed_time=0.0;
double dt=0.002;
double t;
int atomID=0;
int Natoms;
int nmoves=100000;
double beta=0.0;
double T=100; //K
double R=0.0, Rcut=2.5;
double dx_rand, dy_rand, dz_rand, d_max=0.5;
double R2cut=0.0, RIcut=0.0, RI2cut=0.0, RI6cut=0.0, RI12cut=0.0;
double **r=NULL, **r_trial=NULL,  **f=NULL;
double lx=6.8, ly=6.8, lz=6.8;
double Hlx=0.5*lx, Hly=0.5*ly, Hlz=0.5*lz;
double Vir=0.0, Vol=0.0;
double KB=1.0, P=0.0;
double delU=0.0, DU=0.0;
int count_acc=0;
double eps = 0.0;
double acc = 0.0;

// Energy and force calculation variables
double U_init=0.0, U=0.0;
double dx=0.0, dy=0.0, dz=0.0;
double r2=0.0, r6=0.0, Ir6=0.0;
double fcut=0.0, fijx=0.0, fijy=0.0, fijz=0.0, F=0.0;

//-----------------------------------------------------------------//
/* 
 * brief - Write the coordinates of the atom at every timestep
 */
void write_xyz(std::ofstream& simFile, int config)
{
  double Xcom=0.0, Ycom=0.0, Zcom=0.0;
  simFile << "    " << Natoms << std::endl;
  simFile << " i =  " << config << std::endl;
  for(int i=0; i<Natoms; i++) {
    simFile << "Ar" << std::setw(15) << r[i][0] << std::setw(20) << r[i][1] << std::setw(20) << r[i][2] << std::endl;
    Xcom += r[i][0];
    Ycom += r[i][1];
    Zcom += r[i][2];
  }
  //  simFile << "O" << std::setw(15) << Xcom/Natoms << std::setw(20) << Ycom/Natoms << std::setw(20) << Zcom/Natoms << std::endl;
}

//-----------------------------------------------------------------//
/* 
 * brief - Write the stats information
 */
void dump_stats(std::ofstream& enerFile, int config)
{
  enerFile << std::setw(8) << config
	   << std::setw(12) << std::setw(12) << P
	   << std::setw(12) << U  << std::endl;
}

//-----------------------------------------------------------------//
/*
 * brief - get number of atoms
 */
int get_natoms(const char* filename)
{
  int count=0;
  double x,y,z; // few dummy variables
  std::ifstream in_xyz;
  in_xyz.open(filename, std::ifstream::in);

  while( !in_xyz.eof() ){
    in_xyz >> x >> y >> z;
    count++;
  }
  in_xyz.close();
  return (count-1);
}

//-----------------------------------------------------------------//
/*
 * brief - read the coordinates of atoms from a text file
 */
void read_xyz(const int Natoms, const char* filename)
{
  int count=0;
  std::cout << "Loading coordinates.." << std::endl;
  std::ifstream in_xyz;
  in_xyz.open(filename,std::ifstream::in);

  while(count < Natoms){
    in_xyz >> r[count][0] >> r[count][1] >> r[count][2];
    count++;
  }
  in_xyz.close();
}

//-----------------------------------------------------------------//
/*
 * brief - Select an atom randomly [0-255] and perturb it in 
 *         x, y, z directions by dx,dy,dz.
 */
void compute_trial_params(){
  
  // Set the perturbations (-1,1)
  dx_rand = d_max * (2.0*(double(rand())/double(RAND_MAX)) - 1.0);
  dy_rand = d_max * (2.0*(double(rand())/double(RAND_MAX)) - 1.0);
  dz_rand = d_max * (2.0*(double(rand())/double(RAND_MAX)) - 1.0);

  // Randomly select an atom (Just to make sure that the rand no is not same as the last one)
  int temp_atomID=0;
  bool foundatom=false;
  while(!foundatom){
    temp_atomID = rand()%Natoms;
    if(temp_atomID != atomID) {
      foundatom = true;
      atomID = temp_atomID;
    }
  }
  if(atomID > Natoms){
    std::ostringstream msg;
    msg << __LINE__ << " : " << __FILE__ << std::endl
	<< "Error ! Trying to randomly move an atom that doesn't exist. [atomID] :"  << atomID << std::endl;
    throw std::runtime_error( msg.str() );
  }
  
  for(int i=0; i<Natoms; i++){
    if( atomID==i ){
      r_trial[atomID][0] = r[atomID][0] + dx_rand;
      r_trial[atomID][1] = r[atomID][1] + dy_rand;
      r_trial[atomID][2] = r[atomID][2] + dz_rand;
    }
    else{
      r_trial[i][0] = r[i][0];
      r_trial[i][1] = r[i][1];
      r_trial[i][2] = r[i][2];
    }

    // Apply PBC
    if(r_trial[i][0] <  0) r_trial[i][0] = r_trial[i][0] + lx;
    if(r_trial[i][0] > lx) r_trial[i][0] = r_trial[i][0] - lx;
    if(r_trial[i][1] <  0) r_trial[i][1] = r_trial[i][1] + ly;
    if(r_trial[i][1] > ly) r_trial[i][1] = r_trial[i][1] - ly;
    if(r_trial[i][2] <  0) r_trial[i][2] = r_trial[i][2] + lz;
    if(r_trial[i][2] > lz) r_trial[i][2] = r_trial[i][2] - lz;
  }
}

//-----------------------------------------------------------------//

/*
 * brief - initialize the positions by reading the inputs.
 *
 */
void init(int Natoms, const char* filename){
  r     = new double*[Natoms];
  r_trial = new double*[Natoms];
  f     = new double*[Natoms];

  for(int i=0; i<Natoms; i++){
    r[i]     = new double[3];
    r_trial[i] = new double[3];
    f[i]     = new double[3];
  }

  // initialize the positions
  read_xyz(Natoms, filename);
}

//-----------------------------------------------------------------//

/* 
 * brief - Employ Metropolis Monte-Carlo
 */
void apply_metropolis()
{

  if(delU<0){ // Accept the move
    // Update the positions
    for(int i=0; i<Natoms; i++){
      r[i][0] = r_trial[i][0];
      r[i][1] = r_trial[i][1];
      r[i][2] = r_trial[i][2];
    }
    U = compute_U();
  }
  else{
    eps = double(rand())/double(RAND_MAX);
    acc = exp(-beta*(U_post-U_pre));
    if(eps <= acc){ // Accept
      count_acc++;

      for(int i=0; i<Natoms; i++){
	r[i][0] = r_trial[i][0];
	r[i][1] = r_trial[i][1];
	r[i][2] = r_trial[i][2];
      }
      U = compute_U();
    }
    else{ // Reject
      U = compute_U();
    }
  }
}

//-----------------------------------------------------------------//
// void compute_preU()
// {
//   U_pre = 0.0;
//   double dx, dy, dz, dr2, dr6, dr12;
//   for(int i=0; i<Natoms; i++){
//     for(int j=i+1; j<Natoms; j++){
//       dx = r[i][0] - r[j][0];
//       dy = r[i][1] - r[j][1];
//       dz = r[i][2] - r[j][2];
//       dr2 = (dx*dx)+(dy*dy)+(dz*dz);
//       dr6 = dr2*dr2*dr2;
//       dr12= dr6*dr6;

//       U_pre += 4.0*((1.0/dr12) - (1.0/dr6));
//     }
//   }
// }

// void compute_postU()
// {
//   U_post=0.0;
//   double dx, dy, dz, dr2, dr6, dr12;
//   for(int i=0; i<Natoms; i++){
//     for(int j=i+1; j<Natoms; j++){
//       dx = r_trial[i][0] - r_trial[j][0];
//       dy = r_trial[i][1] - r_trial[j][1];
//       dz = r_trial[i][2] - r_trial[j][2];
//       dr2 = (dx*dx)+(dy*dy)+(dz*dz);
//       dr6 = dr2*dr2*dr2;
//       dr12= dr6*dr6;

//       U_post += 4.0*((1.0/dr12) - (1.0/dr6));
//     }
//   }
// }

double compute_U()
{
  U_init = 0.0;
  double dx=0.0, dy=0.0, dz=0.0, dr2=0.0, Idr6=0.0;
  for(int i=0; i<Natoms; i++){
    for(int j=i+1; j<Natoms; j++){
      dx = r[i][0] - r[j][0];
      dy = r[i][1] - r[j][1];
      dz = r[i][2] - r[j][2];
      dr2 = (dx*dx)+(dy*dy)+(dz*dz);
      Idr6 = 1/(dr2*dr2*dr2);

      U_init += 4.0*(Idr6*Idr6 - Idr6);
    }
  }
  return U_init;
}

//-----------------------------------------------------------------//

void compute_delU()
{
  // abhi : There is some problem with delU
  double U_trial=0.0, U_old=0.0;
  double dx_trial,dy_trial,dz_trial,dr2_trial,dr6_trial,dr12_trial;
  double dx,dy,dz,dr2,dr6,dr12;;
  delU = 0.0;

  for(int i=0; i<Natoms; i++){
    if(atomID!=i){
      dx_trial = r_trial[i][0] - r_trial[atomID][0];
      dy_trial = r_trial[i][1] - r_trial[atomID][1];
      dz_trial = r_trial[i][2] - r_trial[atomID][2];
      dr2_trial = (dx_trial*dx_trial)+(dy_trial*dy_trial)+(dz_trial*dz_trial);
      dr6_trial = dr2_trial*dr2_trial*dr2_trial;
      dr12_trial= dr6_trial*dr6_trial;
      
      dx   = r[i][0] - r[atomID][0];
      dy   = r[i][1] - r[atomID][1];
      dz   = r[i][2] - r[atomID][2];
      dr2  = (dx*dx)+(dy*dy)+(dz*dz);
      dr6  = dr2*dr2*dr2;
      dr12 = dr6*dr6;

      U_trial = 4.0*((1.0/dr12_trial) - (1.0/dr6_trial));
      U_old   = 4.0*((1.0/dr12)       - (1.0/dr6));
      delU += (U_trial - U_old);
      // if(delU>50){
      // 	std::cout << "Exiting ....U_trial & U_old : " << U_trial << ", " << U_old << ", " << dr12_trial << ", " << dr6_trial << ", " << dr12 << " , " << dr6 << std::endl;
      // 	exit(1);
      // }
    }
  }
  //std::cout << "Value of delU : " << delU << std::endl;
  //DU += delU;
}

//-----------------------------------------------------------------//

void calc_virial_force()
{
  // Before calculating energy, make sure to rezero them
  dx=0.0; dy=0.0; dz=0.0;
  r2 = 0.0; r6=0.0; Ir6=0.0;
  fcut=0.0; F=0.0; Vir=0.0;
  
  // Force cutoff correction : Scheme 2
  fcut = RIcut*(-48*RI12cut + 24*RI6cut);  // dU/dr at Rcut
  
  for(int i=0; i<Natoms; i++){
    f[i][0] = 0.0;
    f[i][1] = 0.0;
    f[i][2] = 0.0;
  }

  for(int i=0; i<Natoms; i++){
    for(int j=i+1; j<Natoms; j++){
      fijx=0.0; fijy=0.0; fijz=0.0;
      R = 0;
      dx = r[i][0] - r[j][0];
      dy = r[i][1] - r[j][1];
      dz = r[i][2] - r[j][2];

      // Apply PBC: Nearest Image Convention
      if(dx > Hlx) dx=dx-lx;
      if(dx <-Hlx) dx=dx+lx;
      if(dy > Hly) dy=dy-ly;
      if(dy <-Hly) dy=dy+ly;
      if(dz > Hlz) dz=dz-lz;
      if(dz <-Hlz) dz=dz+lz;

      r2 = (dx*dx) + (dy*dy) + (dz*dz);

      if(r2<R2cut){
	r6 = r2*r2*r2;
	Ir6 = 1.0/r6;
	R = std::sqrt(r2);

	F    = 24*(2*Ir6*Ir6 - Ir6);
	fijx = F*(dx/r2) + fcut*(dx/R);
	fijy = F*(dy/r2) + fcut*(dy/R);
	fijz = F*(dz/r2) + fcut*(dz/R);
	
	f[i][0] += fijx;
	f[i][1] += fijy;
	f[i][2] += fijz;
	f[j][0] -= fijx;
	f[j][1] -= fijy;
	f[j][2] -= fijz;

	Vir += (F + (fcut*R));
      } // r2< R2cut
    }
  }
  //std::cout << "F:  " << F << std::endl;
}

//-----------------------------------------------------------------//

void calc_inst_pressure()
{
  P=0.0;
  P = ((Natoms*T)/(Vol*120.92) + (Vir/(3*Vol))) * 42.49; //MPa
  //std::cout << "Value of temp : " << temp << std::endl;
}

//-----------------------------------------------------------------//

int main(int argc, char** argv)
{
  // Load the text file
  const char filename[]="liquid256.txt";

  // Determine the number of atoms
  Natoms = get_natoms(filename);

  // Initialize the positions !
  init(Natoms, filename);

  // intialize the random seed !
  srand (time(NULL));

  // Compute some parameters
  R2cut   = Rcut*Rcut;
  RIcut   = 1/Rcut;
  RI2cut  = 1/R2cut;
  RI6cut  = RI2cut*RI2cut*RI2cut;
  RI12cut = RI6cut*RI6cut;
  Vol     = lx*ly*lz;
  beta    = 1/(KB*T);

  std::ofstream simFile,enerFile;
  simFile.open("LDmj_sim.xyz");
  enerFile.open("LDmj_sim.ener");
  std::cout << "----------------------------------------" << std::endl;
  std::cout<< "     move" << std::setw(15) << std::setw(12) << "P"
	   << std::setw(12) << "U" << std::endl;
  std::cout << "----------------------------------------" << std::endl;

  
  // Loop over trial moves
  for(int k=0; k<nmoves; k++){
    //if(k==0) compute_U();

    //compute_preU();
    
    // Compute trial moves
    compute_trial_params();

    //compute_postU();
    compute_delU();

    // Apply Metropolis MC
    apply_metropolis();

    calc_virial_force();
    calc_inst_pressure();

    // Output the results for every 100 moves
    if(k%100 == 0){
      std::cout << std::setw(8)  << k  << std::setw(12)
      		<< std::setw(25) << P  << std::setw(15)
      		<< std::setw(25) << U  << std::setw(20) << delU << std::endl;
      
      write_xyz(simFile, k);
      dump_stats(enerFile, k);
    }
  }
  
  simFile.close();
  enerFile.close();
  
  return 0;
}
