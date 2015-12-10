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

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <iomanip>
#include <stdexcept>

// Function definitions
void dump_stats(std::ofstream& enerFile, int config);
double compute_U();
void compute_trial_params();
void apply_metropolis();

//-----------------------------------------------------------------//
// Energy and force calculation variables
double U_trial=0.0, U=0.0;
double dx_max=5, dx=0.0;
double x_trial=0.0, x=0.0;

double q=0.0;
double beta=0.01;
//double BETA[5]={0.01, 0.1, 1.0, 10.0, 100.0};
//double DX_MAX[5]={0.01, 0.05, 0.1, 0.5, 5};
double EPSILON[5]={0.01, 0.05, 0.1, 0.005};
double epsilon=0.005;
double pos_end=q+(0.5*epsilon);
double neg_end=q-(0.5*epsilon);
int M=0;
int nmoves=10000000;
int count_acc = 0;

//-----------------------------------------------------------------//
// brief - Write the stats information
void dump_stats(std::ofstream& enerFile, int config)
{
  enerFile << std::setw(15) << config
	   << std::setw(20) << U  << std::setw(20) << x 
	   << std::endl;
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
    if(neg_end<x && x<pos_end) M++;
  }
  else{
    B = double(rand())/double(RAND_MAX);
    acc = exp(-beta*(U_trial - U));
    if(B <= acc){ // Accept
      x = x_trial;
      U = U_trial;
      count_acc++;
      if(neg_end<x && x<pos_end) M++;
    }
    else{ // Reject
      x = x;
      U = U;
    }
  }
}

//-----------------------------------------------------------------//

void compute_trial_params()
{
  dx=0.0;
  // Generate trial move
  double start=-2.0;
  double x_start = start + ((double(rand())/double(RAND_MAX)) * (pos_end-start));
  if(start > x_start || x_start > pos_end){
    std::ostringstream msg;
    msg << __LINE__ << " : " << __FILE__ << std::endl
	<< "Error ! Check x_start !!!" << std::endl;
    throw std::runtime_error(msg.str());
  }
  dx = dx_max * (2.0*(double(rand())/double(RAND_MAX)) - 1.0);

  // Generate trial move
  x_trial = x_start + dx;

  // Generate trial potential
  if(count_acc == 0) x = x_start;
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

  // Loop over trial moves
  for(int j=0; j<4; j++){
    //std::ofstream enerFile;
    //std::stringstream filename;
    //filename << "K_HTST_NVT_MC-dx_max-" << DX_MAX[j] << ".txt";
    //enerFile.open(filename.str().c_str());
    //beta=BETA[j];
    //dx_max = DX_MAX[j];
    epsilon=EPSILON[j];
    U=0.0; U_trial=0.0;
    x=0.0; M=0; count_acc=0;

    for(int k=0; k<nmoves; k++){
      // Compute trial moves
      compute_trial_params();
      
      // Apply Metropolis MC
      apply_metropolis();
      
      // Output the results for every 100 moves
      //if(k%100==0) dump_stats(enerFile, k);
    }

    //enerFile.close();
    std::cout << "% Acceptance of trail moves : " << (double(count_acc)/double(nmoves))*100.0 
	      << " , M : " << M << std::endl;
    double temp = 2.0/(M_PI*beta);
    double prob = double(M)/(double(nmoves)*epsilon);
    std::cout << "Value of k(Mc, " << beta << ") : " << 0.5*std::sqrt(temp)*prob << std::endl;
  } // beta

  return 0;
}
