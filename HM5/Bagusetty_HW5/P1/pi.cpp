#include<iostream>
#include<math.h>
#include<fstream>
#include<stdlib.h>
#include<iomanip>
#include<ctime>

long double pi =0.0;
long double eps=0.0;  // estimated uncertainity
double xlow=-1.0, xhigh=1.0, ylow=-1.0, yhigh=1.0;

// ----------------------------------------------------------------------------

int main(){

  int Ndarts=10;   // maximum value of random numbers for producing HITs.
  double x=0.0, y=0.0;      // Coordinates
  double z;
  double xrange = xhigh - xlow;
  double yrange = yhigh - ylow;

  srand(time(NULL));
  double sum=0.0;
  
  for(int i=0; i<Ndarts; i++){
    // produce the values between -1.0 : 1.0
    x = 2.0*(double(rand())/double(RAND_MAX)) - 1;
    y = 2.0*(double(rand())/double(RAND_MAX)) - 1;
    z = x*x + y*y;
    
    if(z<=1) sum = sum + 1.0;
    else     sum = sum + 0.0;

    pi = double(sum/i) * xrange * yrange;
  }

  // Compute the value of PI and it's deviation
  eps = fabs(pi - M_PI);              // estimated uncertainity

  std::cout << std::setprecision(10) << "  " << pi << " , " << eps << std::endl;
  
}


