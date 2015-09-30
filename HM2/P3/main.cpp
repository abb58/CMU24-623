/*
 * Molecular Dynamics simulation of lennard-jones nano-
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
 * Loop unrolling is implemented for the coordinates so
 * hardware thread gets appropriate work.
 *
 * (c) 2015-16
 */

#include <iostream>
#include <fstream>
#include <cmath>
#include "driver.h"
#include "fileio.cpp"


//-----------------------------------------------------------------//

/*
 * brief - initialize the positions by reading the inputs.
 *
 */
void init(int Natoms, const char* filename){

  // Memory Management
  m     = new double[Natoms];
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

  // intialize the mass, velocities, forces
  for(int i=0; i<Natoms; i++){
    m[i] = 1.0;
    v[i][0] = 0.0; v[i][1] = 0.0; v[i][2] = 0.0;
    f[i][0] = 0.0; f[i][1] = 0.0; f[i][2] = 0.0;
  }
}

//-----------------------------------------------------------------//

void vv_scheme()
{
  double tfact;
  
  // TODO : Perform openMP.
  // reinitialize all the temp variables to zero.
  for(int i=0; i<Natoms; i++){
    r_old[i][0] = r[i][0];
    r_old[i][1] = r[i][1];
    r_old[i][2] = r[i][2];

    v_old[i][0] = v[i][0];
    v_old[i][1] = v[i][1];
    v_old[i][2] = v[i][2];

    f_old[i][0] = f[i][0];
    f_old[i][1] = f[i][1];
    f_old[i][2] = f[i][2];
  }

  for(int i=0; i<Natoms; i++){
    tfact = dt/(2.0*m[i]);
    // Step 1. v(t+0.5dt)
    v[i][0] = v_old[i][0] + f_old[i][0]*tfact;
    v[i][1] = v_old[i][1] + f_old[i][1]*tfact;
    v[i][2] = v_old[i][2] + f_old[i][2]*tfact;

    v_old[i][0] = v[i][0];
    v_old[i][1] = v[i][1];
    v_old[i][2] = v[i][2];

    // Step 2. r(t+dt)
    r[i][0] = r_old[i][0] + v_old[i][0]*dt;
    r[i][1] = r_old[i][1] + v_old[i][1]*dt;
    r[i][2] = r_old[i][2] + v_old[i][2]*dt;
  }

  calc_pairenergy();
  calc_force();

  for(int i=0; i<Natoms; i++){
    // Step 3. v(t+dt)
    v[i][0] = v_old[i][0] + f[i][0]*tfact;
    v[i][1] = v_old[i][1] + f[i][1]*tfact;
    v[i][2] = v_old[i][2] + f[i][2]*tfact;
  }
}

//-----------------------------------------------------------------//

void calc_pairenergy()
{
  // Before calculating energy, make sure to rezero them
  U = 0.0;
  dx=0.0; dy=0.0; dz=0.0;
  r2 = 0.0; r6=0.0; Ir6=0.0;
  
  for(int i=0; i<Natoms; i++){
    for(int j=i+1; j<Natoms; j++){
      dx = r[i][0] - r[j][0];
      dy = r[i][1] - r[j][1];
      dz = r[i][2] - r[j][2];
      r2 = (dx*dx) + (dy*dy) + (dz*dz);
      r6 = r2*r2*r2;
      Ir6 = 1/r6;

      // Pair Energy
      U += 4*(Ir6*Ir6 - Ir6);
    }
  }
}

//-----------------------------------------------------------------//

void calc_force()
{
  // Before calculating energy, make sure to rezero them
  dx=0.0; dy=0.0; dz=0.0;
  r2 = 0.0; r6=0.0; Ir6=0.0;
  F = 0.0;
  
  for(int i=0; i<Natoms; i++){
    f[i][0] = 0.0;
    f[i][1] = 0.0;
    f[i][2] = 0.0;
  }

  for(int i=0; i<Natoms; i++){
    for(int j=i+1; j<Natoms; j++){
      dx = r[i][0] - r[j][0];
      dy = r[i][1] - r[j][1];
      dz = r[i][2] - r[j][2];
      r2 = (dx*dx) + (dy*dy) + (dz*dz);
      r6 = r2*r2*r2;
      Ir6 = 1.0/r6;

      // LJ Force
      // This procedure of manipulating r avoids use of sqrt()
      F        = 24*(2*Ir6*Ir6 - Ir6);
      f[i][0] += F * (dx/r2);
      f[i][1] += F * (dy/r2);
      f[i][2] += F * (dz/r2);

      f[j][0] -= F * (dx/r2);
      f[j][1] -= F * (dy/r2);
      f[j][2] -= F * (dz/r2);
    }
  }
}

//-----------------------------------------------------------------//

void calc_kenergy()
{
  KE = 0.0;
  for(int i=0; i<Natoms; i++) KE += (v[i][0]*v[i][0]) + (v[i][1]*v[i][1]) + (v[i][2]*v[i][2]);
  KE *= 0.5;
}

//-----------------------------------------------------------------//

void calc_momentum()
{
  px=0.0; py=0.0; pz=0.0;
  // Assumption : Mass of atoms = 1.0
  for(int i=0; i<Natoms; i++){
    px += m[i]*v[i][0];
    py += m[i]*v[i][1];
    pz += m[i]*v[i][2];
  }
}

//-----------------------------------------------------------------//

int main(int argc, char** argv)
{
  const char filename[]="10.txt";

  // Determine the number of atoms
  Natoms = get_natoms(filename);

  // Initialize the positions !
  init(Natoms, filename);

  // Loop over timesteps
  std::ofstream simFile,enerFile;
  simFile.open("LDmj_sim.xyz");
  enerFile.open("LDmj_sim.ener");
  enerFile << "Time step" << std::setw(15) << "time"
	   << std::setw(15) << "PE" << std::setw(15) << "KE"
	   << std::setw(15) << "TE" << std::setw(15) << "Px"
	   << std::setw(15) << "Py" << std::setw(15)
	   << "Pz" << std::setw(15) << std::endl;

  for(int k=0; k<=1000; k++) {
    elapsed_time = dt*double(k);

    // Calculate pair-energy and forces
    if(k==0) calc_force();

    calc_kenergy();
    TE=U+KE;

    std::cout << std::setw(8)  << k  << std::setw(15) << elapsed_time
	      << std::setw(15) << U  << std::setw(15) << KE
	      << std::setw(15) << TE << std::setw(15) << px
	      << std::setw(15) << py << std::setw(15) << pz
	      << std::setw(15) << std::endl;

    calc_momentum();
    write_xyz(simFile, k);
    dump_stats(enerFile, k);

    // parameters are computed for (i+1)
    vv_scheme();
  }
  
  simFile.close();
  enerFile.close();
  
  return 0;
}
