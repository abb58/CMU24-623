/*
 * Molecular Dynamics simulatio of lennard-jones nano
 * particle.
 * 
 * Abhishek Bagusetty
 * 
 * Writte for CMU24-623 Molecular Simulation of Materials,
 * Fall 2015-16
 *
 * Building instructions :
 *    - make
 *    - make clean
 * 
 * External Dependencies :
 * GNU Scientific Library is required.
 * 
 * (c) 2015-16
 */

#include <iostream>
#include <cmath>
#include "fileio.cpp"

//------ Function Declarations ------
void write_xyz();
void read_xyz();
void print_energy();
void init();
double tot_energy();
double pot_energy();
double kin_energy();
void integrate();

//------ Function Definitions ------
/*
 * brief - write enrgy file
 */
void print_energy()
{
}



int main(int argc, char** argv)
{
  double **r, **v, **f;
  


  return 0;
}
