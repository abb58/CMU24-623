// start.cpp
//
// This program performs a series of mathematical operations on a set of integers read in from a file.
// 
//*****************
#include <iostream> // this is a standard C++ header, and should be included in all programs
#include <fstream>
#include <iomanip>
#include <math.h>
#include <time.h>
using namespace std;
//*****************

// declaration of functions/subroutines
int i_polynomial(int i);	// int - a polynomial to be evaluated
double d_polynomial(int i);	// double - a polynomial to be evaluated


// declaration of input and output stream
ifstream data_in("5input.txt"); // input data
ofstream data_out("5output.txt"); // output data


// the main part of the program - should always be type "int"

int main() {
	
  // variable declarations
  int i,j,k; // loop counters
  const int upper = 11; // upper limit on calculation loop, corresponds to number of entries in the input file
    					  // declared as a "const int" so that arrays can be defined with this size
  int numbers[upper]; // array to store numbers, in C++, the array index starts at 0!
    
  bool sq;          // boolean that indicates if the number is a perfect square
  int pol;          // value of polynomial function
  int quotient, divisor, remainder;
  char square[200]; // strings for output
  int sqrtN;
  bool purge;       // condition to purge the for loop 
  
  // read in data from file
  for(i=0;i<upper;i++) data_in>>numbers[i];
  // numbers[0] = 11;
  
  // calculation loop
  for(i=0;i<upper;i++){ // the loop will run over the elements in numbers[]    

    // first, determine if the number is a perfect square
    sprintf(square," is not a perfect square"); //default
    if(numbers[i]<0) sprintf(square," is negative so that the concept of a perfect square is undefined for real numbers");
    else {
      quotient = numbers[i];
      divisor = 2;
      sq = false;
      while(sq == false && quotient > divisor) {
	quotient  = numbers[i]/divisor;
	remainder = numbers[i]%divisor;
	if (quotient == divisor && remainder == 0) {
	  sq = true;
	  sprintf(square," is a perfect square");
	  sqrtN = sqrt(numbers[i]);
	}
	else divisor = divisor + 1;
      }
    }

    // Check if the integer is a PRIME
    // Primality Algorithm
    for(j=2; j<numbers[i]; j+=2) { // i - divisor

      // Divisible by 2 ?
      if( j==2 ){
	remainder = numbers[i]%j;
	if(remainder==0) {
	  data_out << "The integer " << numbers[i] << " is COMPOSITE" << endl;
	  cout << "The integer " << numbers[i] << " is COMPOSITE" << endl;
	  break;
	}
	else{
	  j+=1;
	  remainder = numbers[i]%j;
	  if(remainder==0) {
	    data_out << "The integer " << numbers[i] << " is COMPOSITE" << endl;
	    cout << "The integer " << numbers[i] << " is COMPOSITE" << endl;
	    break;
	  }
	}
      }

      if(sq && j<=sqrtN){ // perfect square
	remainder = numbers[i]%j;
	if(remainder == 0){
	  data_out << "The integer " << numbers[i] << " is COMPOSITE" << endl;
	  cout << "The integer " << numbers[i] << " is COMPOSITE" << endl;
	  break;
	}
	else if(j == sqrtN){
	  data_out << "The integer " << numbers[i] << " is PRIME" << endl;
	  cout << "The integer " << numbers[i] << " is PRIME" << endl;
	  break;
	}
	else{}
      }
      else{ // Non-Perfect Square ?
	remainder = numbers[i]%j;
	if(remainder==0) {
	  data_out << "The integer " << numbers[i] << " is COMPOSITE" << endl;
	  cout << "The integer " << numbers[i] << " is COMPOSITE" << endl;
	  break;
	}
	else if(j == numbers[i]-2){
	  data_out << "The integer " << numbers[i] << " is PRIME" << endl;
	  cout << "The integer " << numbers[i] << " is PRIME" << endl;
	  break;
	}
	else{};
      }

    } // PRIME check loop
    
    // second, evaluate the polynomial function
    pol = i_polynomial(numbers[i]);
    pol = d_polynomial(numbers[i]);
    
    //output the data
    //cout<<"The number "<<numbers[i]<<square<<" and returns a value of "<<pol<<" when inserted into the function f"<<endl;  
  } 
}	

/*
 * \brief : method to evaluate a polynomial and return the output in int datatype
 */
 ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int i_polynomial(int i) {
  int x2 = i*i;
  int x4 = x2 * x2;
  int x6 = x4 * x2;
  int f;
  f = 2*x6 - 3*x4 + 4*x2 - 3;
  return f;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*
 * \brief : method to evaluate a polynomial and return the output in double datatype
 */
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double d_polynomial(int i) {
  double x2 = i*i;
  double x4 = x2 * x2;
  double x6 = x4 * x2;
  double f;
  f = 2.0*x6 - 3.0*x4 + 4.0*x2 - 3.0;
  return f;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

