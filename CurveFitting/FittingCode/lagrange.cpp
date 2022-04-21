// Copyright 2022 Richard Brower brower@bu.edu
#include <iostream>
#include <fstream>
#include <chrono>
#include <cmath>
using namespace std;

/*****
see\href{https://en.wikipedia.org/wiki/Overfitting}{https://en.wikipedia.org/wiki/Overfitting}

\href{https://en.wikipedia.org/wiki/Lagrange\_polynomial}{https://en.wikipedia.org/wiki/Lagrange\_polynomial}

******/

double fExact(double x, double *xP, double *yP, int Npoint)
{
  double f, lagrange = 0.0;
  f = 0;
  for(int i = 0; i <Npoint ; i++)
    {
      lagrange = 1.0;	  
      for(int j = 0; j< Npoint ;  j++)
	{
	  if(j != i) lagrange *= ( x - xP[j] )/( xP[i] - xP[j] );
	}
      f += yP[i]*lagrange;
    }  
  return f;
}

int main()
{
  int Npoint = 20;
  double x;
  double xP[Npoint];
  double yP[Npoint];
  double ranf;  // in range of [0,1]
  // Make data    Fdata[NPoint][2]


  for(int i = 0; i <Npoint ; i++)
    {
      ranf = (double)rand()/(double)RAND_MAX; 
      xP[i] = 10.0 * (i+  0.2* ranf)/(double)Npoint ;
      ranf = (double)rand()/(double)RAND_MAX; 
      yP[i] = (1.0 + 2.0 *xP[i] + 3.0 * xP[i]*xP[i]) * (1.0 +  0.5 * ( ranf - 0.5));					  
        cout << "   " <<  xP[i] << "   "  <<  yP[i] <<  endl;
    };

  // double f = fExact(2.0, xP, yP, Npoint);
  //cout <<"  Test fExactfExact(x, xP, yP, Npoint) "<< "   "<<  f  << endl;
  
  cout << endl  << "# The Fitted Form " << endl;

  #if 1
  int refinement = 10;
  for(int i = 0; i < refinement*Npoint ; i++)
    {
      x =  xP[0] + i*(xP[Npoint-1] -  xP[0])/(refinement*(double)Npoint);
      double ff = fExact(x, xP, yP, Npoint);
      //  cout << " i =    " << i << " x =   " << x << "   " <<  ff << endl;
       cout << "  " << x << "   " <<  ff << endl;
    }
 #endif 
  return 0;
}
