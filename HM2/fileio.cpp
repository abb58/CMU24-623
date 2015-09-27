#include <fstream>
#include <iostream>

/* 
 * brief - Write the coordinates of the atom at every timestep
 */
void write_xyz(std::ofstream& simFile, int config)
{
  double Xcom=0.0, Ycom=0.0, Zcom=0.0;
  simFile << "    " << Natoms+1 << std::endl;
  simFile << " i =  " << config << std::endl;
  for(int i=0; i<Natoms; i++) {
    simFile << "H       " << r[i][0] << "      " << r[i][1] << "      " << r[i][2] << std::endl;
    Xcom += r[i][0];
    Ycom += r[i][1];
    Zcom += r[i][2];
  }
  simFile << "O       " << Xcom/Natoms << "      " << Ycom/Natoms << "      " << Zcom/Natoms << std::endl;
}


/* 
 * brief - Write the stats information
 */
void dump_stats(std::ofstream& enerFile, int config)
{
  enerFile << "   "  << config << "      " << elapsed_time << "      " << U << "       " << KE
	   << "       " << TE << "         " << px << "       " << py << "       " << py << std::endl;
}


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
