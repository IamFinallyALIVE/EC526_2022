/*************************
This code calculate the N-th root  of A  which is  x = A^(1/N)  
e.g. in C  double pow(double x, double y);  x = pow(A,1.0/N) 



The intialization sets A and N and the search range is
automatically  x = A^(1/N) \in [0, A] for A > 1 and  N = 2,3, 4...

Note we are not cheating because the rootines only use
integer N   double pow(double x,  N);  which could be done with iteration 
to make this clear with log_2(N) multiplication. 

However the pow is only double  with 17 decimal digits for binary64,
to go to long double   with 36 decimal digits for binary128
we would have write our own pow

 ************************/
#include <iostream>
#include <iomanip> 
#include <cmath>
#include <float.h> // For long doubles.
#include <time.h> // second in the Unix Era

using std::cout;
using std::cin;
using std::endl;
using std::ios;

using std::abs;

// This function performs one iteration of Newton's method
// and returns a new guess (x - f(x)/f'(x) -> x_new).

double newton(double x, double A, double N);

// This function performs one iteration of bisection
// and updates either "min" or "max" (note how they are both
// passed by reference), and returns the current "midpoint".
// Again, you need to hard-code the numerical function. Bisection
// does not require a derivative. 
double bisection(double A, double N, double  & min, double  & max);

int main()
{
  // Declare variables to hold the current guess
  // and relative error. 
  double  x = 0.0, fractional_error = 0.5;
  
  // Declare a variable to hold "A". 
  double  A;

  // Declare a root "N" to find
   double N;

  // Declare a tolerance
  double  log10_tol;
   double  tol = 0.0;
  
  // Declare a counter.
  int count; 
  
  // Print the number of digits a double  can hold.
  cout << "Number of digits accuracy in double  " << LDBL_MANT_DIG << endl;
  
  // Fix the output format.
  cout.setf(ios::fixed,ios::floatfield); 


  cout.precision(40);
  
  // Describe the problem and prompt the user:
  cout << " Compute the Nth root by Newton's Method and the bisection method to a tolerance " << tol << "." << endl;
  cout <<"Give a number A: ";
  cin >> A;
  cout << "Give a number N: ";
  cin >> N;
  cout << "Give a log10 tolerance (i.e., tolerance will be 10^{-[number]}): ";
  cin >> log10_tol;
  tol = pow(10.0,-log10_tol);

  // Choose an initial guess for Newton's method: in this case,
  // A/N. Set the output precision as well.
  x = A/N;
  count = 0;
  cout.precision(40);
  do 
  {
    count++;
    x = newton(x, A, N);
    fractional_error = 0.5*abs(pow(x,N)/A-1.0);
    cout << x << "\t" << fractional_error << endl;
  }
  while(fractional_error > tol && count < 100000);

  cout.precision(40);
  cout << "After " << count << " iterations, Newton's method " << endl;
  cout << "gave = "  << x  <<  " vs cmath = "  << pow(A,1.0/N) << endl; 
  cout.precision(40);
  cout << " relative error   = " << (x - pow(A,1.0/N))/(pow(A,1.0/N)) << endl;
  
  /*   Compare with bisection */
  cout<< " Bisection Method starting x = A/N  with min = 0 and max = A \n";
  
  double  min, max;
  min = 0.0;
  max = A;
  count = 0;
  
  do 
    {
      count++;
      x =  bisection( A, N,  min, max);
      fractional_error = 0.5*abs(pow(x,N)/A-1);
      cout.precision(40);
      cout << x   << "\t" <<  fractional_error << endl;
    }
  while(fractional_error > tol);
  
  
  cout.precision(40);
  cout << "Bisection's value  in " << count << " iterations " << endl;
  cout << "gave = "  << x  <<  " vs cmath = "  << pow(A, 1.0/N) << endl;
  cout.precision(40);
  cout << " error   = " << x - pow(A, 1.0/N) << endl;
  
  //  Make file for  table for graphing using gnuplot
  
  /*** Simplest  4 steps to write to a file****/
  // FILE * fptr;
  // fptr = fopen("rootData.txt","w");
  // INSERT STUFF fprint(fptr,"Blaa Blaa ...."
  // close(fptr);
  /**  Below adding infomation to file such as  trucated UNIX  time stamp **/  

  FILE * fptr;
  char outfileName[50];
  time_t seconds;   
  seconds = time(NULL);
  sprintf(outfileName, "rootData_%ld.dat",seconds%100000);
  fptr = fopen(outfileName,"w");
  int Nint = N;
  fprintf(fptr, "#Data File for N =  %d-th root of A =  %f \n", Nint , A);
  fprintf(fptr, "#tolerance   Bisection Count  Newton Count \n");
  
   tol = 1.0;
  int  Bicount;
  int Newcount;
  
  for(int iter = 0; iter < 15; iter++)
    {
      tol *= 0.1;
      Bicount =0;
      min = 0.0;
      max = A;
      x = A/N;
      do 
	{
	  Bicount++;
	  x =  bisection( A, N,  min, max);
	  fractional_error = 0.5*abs(pow(x,N)/A-1.0);
	}
      while(fractional_error > tol&& Bicount < 100000);
      
      Newcount = 0;
      min = 0.0;
      max = A;
      x = A/N;
      do 
	{
	  Newcount++;
	  x = newton(x, A, N);
	  fractional_error = 0.5*abs(pow(x,N)/A-1.0); 
	}
      while(fractional_error > tol && Newcount < 100000);
      
      fprintf(fptr, "  %25.20e         %10d       %10d  \n", 1.0/tol,Bicount, Newcount);
    }
  
  fclose(fptr);


  return  0;
}

// This routine is currently hard coded for the function
// f(x) = x^N-A. or x = x - f(x)/f'(x) = x - (x^N - A)/N x^(N-1)
// or x = x - f(x)/f'(x) = x - (x^N - A)/N x^(N-1) = x*(1 - 1.0/N) +  A/(N*x^(N-1))
//
double  newton(double  x, double  A, double  N)
{
  return x*(1 - 1.0/N) + A/(N*pow(x,N-1.0));
}

// This routine is currently hard coded for the function
// f(x) - x^2 - A
double  bisection(double  A, double  N, double  & min, double  & max)
{
  double  x  = (min + max)/2.0;
  if(pow(x,N)-A < 0.0)
    min = x;
  else
    max = x;

  return x;
}
