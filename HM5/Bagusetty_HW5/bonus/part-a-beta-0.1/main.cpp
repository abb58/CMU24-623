/*
 * Molecular Dynamics simulation of lennard-jones nano-
 * particles.
 *
 * Abhishek Bagusetty
 *
 * Written for CMU24-623 Molecular Simulation of Materials,
 * Nose-Hoover Thermostat
 * Fall 2015-16
 *
 * Building instructions :
 *    - make
 *    - make clean
 * Loop unrolling is implemented for the coordinates so
 * hardware thread gets appropriate work.
 *
 * (c) 2015-16
 */

#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <iomanip>

// Function definitions
void dump_stats(std::ofstream& enerFile, int config);
void compute_trial_parameters();
void apply_metropolis();
void apply_kawasaki();



int n_moves = 100000;

// Energy and force calculation variables
double U_trial=0.0, U=0.0;
double dx_max=5, dx=0.0;
double x_trial = 0.0, x = 0.0;

double eps = 0.0;
double acc = 0.0;
double beta = 0.1;

int count_acc = 0;

//-----------------------------------------------------------------//

/* 
 * brief - Write the stats information
 */
void dump_stats(std::ofstream& enerFile, int config)
{
  enerFile << std::setw(15) << config
	   << std::setw(15) << U  << std::setw(15) << x << std::setw(15) << x*x 
	   << std::endl;
}

//-----------------------------------------------------------------//

/* 
 * brief - Employ Metropolis Monte-Carlo
 */
void apply_metropolis()
{
  if(U_trial < U){ // Accept the move
    x = x_trial;
    U = U_trial;
  }
  else{
    eps = double(rand())/double(RAND_MAX);
    acc = exp(-beta*(U_trial - U));

    if(eps <= acc){ // Accept
      count_acc = count_acc + 1;
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

/* 
 * brief - Employ Kawasaki Monte-Carlo
 */
void apply_kawasaki()
{
  double dE = U_trial - U;
  double A  = exp(-beta*dE*0.5);
  double B  = exp( beta*dE*0.5);

  eps = double(rand())/double(RAND_MAX);
  acc = A/(A+B);

  if(eps <= acc){ // Accept
    count_acc = count_acc + 1;
    x = x_trial;
    U = U_trial;
  }
  else{ // Reject
    x = x;
    U = U;
  }
}

//-----------------------------------------------------------------//

void compute_trial_parameters()
{
  // Generate trial move
  x_trial = x + dx;

  // Generate trial potential
  U       = 0.5*x*x;
  U_trial = 0.5*x_trial*x_trial;
}

//-----------------------------------------------------------------//

int main(int argc, char** argv)
{
  // Set the random number for the positions
  srand (time(NULL));

  // Loop over #of trial moves
  std::ofstream enerFile;
  enerFile.open("LDmj_sim_dx-5.txt");

  for(int k=0; k<=n_moves; k++) {

    // Set the delta for the position 
    dx = dx_max * (2.0*(double(rand())/double(RAND_MAX)) - 1.0);

    // Compute U_trial
    compute_trial_parameters();

    // Metropolis Algorithm
    //apply_metropolis();

    // kawasaki Algorithm
    apply_kawasaki();

    // Write the output to a file
    if(k%100 == 0) {
      std::cout << std::setw(8)  << k  
		<< std::setw(12) << U   << std::setw(12) << x
		<< std::setw(15) << x*x << std::endl;
      dump_stats(enerFile, k);
    }
  }
  enerFile.close();
  std::cout << "Acceptance percentage : " << (double(count_acc)/double(n_moves))*100 << std::endl;
  return 0;
}
