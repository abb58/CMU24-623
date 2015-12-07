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
void compute_trial_params();
void apply_metropolis();

// Energy and force calculation variables
double U_trial=0.0, U=0.0;
double dx_max=5, dx=0.0;
double x_trial=0.0, x=0.0;

double q=-1.0;
double beta = 0.1;
double epsilon = 0.1;
int count_acc = 0;

//-----------------------------------------------------------------//
/* 
 * brief - Write the stats information
 */
void dump_stats(std::ofstream& enerFile, int config)
{
  enerFile << std::setw(8) << config
	   << std::setw(20) << std::setw(20) << P
	   << std::setw(20) << U  << std::endl;
}

//-----------------------------------------------------------------//
/* 
 * brief - Employ Metropolis Monte-Carlo
 */
void apply_metropolis()
{
  double B = 0.0, acc = 0.0;

  if(U_trial < U){ // Accept the move
    x = x_trial;
    U = U_trial;
    count_acc++;
  }
  else{
    B = double(rand())/double(RAND_MAX);
    acc = exp(-beta*(U_trial - U));
    
    if(B <= acc){ // Accept
      count_acc++;
      x = x_trial;
      U = U_trial;
    }
    else{ // Reject
      x = x;
      U = U;
    }
  }
}

//-----------------------------------------------------------------//

void compute_trial_parameters()
{
  // Generate trial move
  double start=-2.0, end=q+(0.5*epsilon)
  double x_start = (double(rand())/double(RAND_MAX));
  dx = dx_max * (2.0*(double(rand())/double(RAND_MAX)) - 1.0);
  x_trial = x_start + dx;

  // Generate trial potential
  double x2=x*x;
  double x2trial=x_trial*x_trial;
  U       = (x2*x2) - (2.0*x2) + 1.0;
  U_trial = (x2trial*x2trial) - (2.0*x2trial) + 1.0;
}

//-----------------------------------------------------------------//

int main(int argc, char** argv)
{
  // intialize the random seed !
  srand (time(NULL));

  std::ofstream enerFile;
  enerFile.open("LDmj_sim.ener");
  
  // Loop over trial moves
  for(int k=0; k<nmoves; k++){

    // Compute trial moves
    compute_trial_params();

    // compute U (post move)
    U_new = compute_U();

    // Apply Metropolis MC
    apply_metropolis();

    // Output the results for every 100 moves
    if(k%100 == 0){
      std::cout << std::setw(8)  << k  << std::setw(12)
      		<< std::setw(25) << U  << std::setw(15)
      		<< std::setw(25) << x  << std::endl;
      dump_stats(enerFile, k);
    }
  }
  
  enerFile.close();
  std::cout << "% Acceptance of trail moves : " << (double(count_acc)/double(nmoves))*100.0 << std::endl;
  return 0;
}
