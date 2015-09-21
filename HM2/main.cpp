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
#include "driver.h"
#include "fileio.cpp"


//------ Function Definitions ------
/*
 * brief - write information on time step , time, total energy,
 *         potential energy, kinetic energy, elapsed time.
 */
void print_data()
{
}

/*
 * brief - initialize the positions by reading the inputs.
 *
 */
void init(int Natoms, const char* filename){

  // Memory Management
  r     = new double*[Natoms];
  v     = new double*[Natoms];
  f     = new double*[Natoms];
  r_old = new double*[Natoms];
  v_old = new double*[Natoms];
  f_old = new double*[Natoms];
  for(int i=0; i<Natoms; i++){
    r[i]     = new double[3];
    v[i]     = new double[3];
    f[i]     = new double[3];
    r_old[i] = new double[3];
    v_old[i] = new double[3];
    f_old[i] = new double[3];
  }
  
  // initialize the positions
  read_xyz(Natoms, filename);

  // intialize the velocities to zero
  for(int i=0; i<Natoms; i++){
    v[i][0] = 0.0; v[i][1] = 0.0; v[i][2] = 0.0;
    f[i][0] = 0.0; f[i][1] = 0.0; f[i][2] = 0.0;
  }
}

void vv_scheme()
{
  double tfact = dt/(2*m);
  // TODO : Do something for the direction !
  // TODO : Perform openMP.
  for(int i=0; i<Natoms; i++){

    // Step 1. v(t+0.5dt)
    v[][] = v_old[][] + f_old[i][]*tfact;;
    v_old[][] = v[][];
    
    // Step 2. r(t+dt)
    r[][] = r_old[][] + v_old[][]*dt;
    
    // Step 3. v(t+dt)
    v[][] = v_old[][] + f[][]*tfact;
  }

}


int main(int argc, char** argv)
{
  const char filename[]="10.txt";

  // Determine the number of atoms
  Natoms = get_natoms(filename);
  std::cout << "1. natoms :" << Natoms << std::endl;

  // Initialize the positions !
  init(Natoms, filename);

  std::cout << "2. natoms :" << Natoms << std::endl;
  return 0;
}
