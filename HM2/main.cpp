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



//-----------------------------------------------------------------//

/*
 * brief - write information on time step , time, total energy,
 *         potential energy, kinetic energy, elapsed time.
 */
void dump_stats()
{
}

//-----------------------------------------------------------------//

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

//-----------------------------------------------------------------//

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

//-----------------------------------------------------------------//

double calc_energy_force()
{

  // Before calculating energy and forces again make sure to rezero them
  u = 0.0;
  for(int i=0; i<Natoms; i++){
    f[i][0] = 0.0; f[i][1] = 0.0; f[i][2] = 0.0;
  }

  
  for(int i=0; i<Natoms; i++){
    for(int j=i+1; j<Natoms; j++){
      dx = r[i][0] - r[j][0];
      dy = r[i][1] - r[j][1];
      dz = r[i][2] - r[j][2];
      r2 = (dx*dx) + (dy*dy) + (dz*dz);
      r6 = r2*r2*r2;
      Ir6 = 1/r6;
      
      // Pair Energy
      u += 4*(Ir6*Ir6 - Ir6);
      
      // LJ Force
      // This procedure of manipulating r avoids use of sqrt()
      F  = 24*(2*Ir6*Ir6 - Ir6);
      f[i][0] += F * (dx/r2);
      f[i][1] += F * (dy/r2);
      f[i][2] += F * (dy/r2);

      f[j][0] -= F * (dx/r2);
      f[j][0] -= F * (dx/r2);
      f[j][0] -= F * (dx/r2);
    }
  }
  return u;
}

//-----------------------------------------------------------------//

double calc_kenergy()
{
  for(int i=0; i<Natoms; i++){
    ke += (v[i][0]*v[i][0]) + (v[i][1]*v[i][1]) + (v[i][2]*v[i][2]);
  }
  ke *= 0.5;
}

//-----------------------------------------------------------------//

int main(int argc, char** argv)
{
  const char filename[]="10.txt";

  // Determine the number of atoms
  Natoms = get_natoms(filename);

  // Initialize the positions !
  init(Natoms, filename);

  // Calculate energy and forces
  U=calc_energy_force();
  KE=calc_kenergy();
  TE=U+KE;
  
  return 0;
}
