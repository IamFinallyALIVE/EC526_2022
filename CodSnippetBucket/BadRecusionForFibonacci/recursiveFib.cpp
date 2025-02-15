/************

Google 

https://en.wikipedia.org/wiki/Fibonacci_number

https://stackoverflow.com/questions/34556986/wrong-output-in-recursive-fibonacci-number-program

Call recusive max of 50 with 40,730,022,147 calls 
unless you have lot of time to waste

https://en.wikipedia.org/wiki/C_data_types

unsigned long long int [0,18,446,744,073,709,551,615]
**************/

#include <iostream> // need to cin and cout
#include <stdlib.h>

using namespace std;



unsigned long long int r = 0;

unsigned long int fib(int n) {
    printf("k: %llu fib n: %d", r++, n);
    if (n==0 || n==1) {
        printf("\n");
        return 1;
    } else { 
        printf(" +\n");
        return fib(n-1) + fib(n-2); 
    }
}

unsigned long int fibClean (int n) {
     r++;
    if (n==0 || n==1) return 1;
    else return fibClean(n-1) + fibClean(n-2);
}


unsigned long int fibIterate(int n) {
  double fib[n+1];
  fib[0] = 1;
  fib[1] = 1;
  int iter;
  for(iter = 2; iter < n+1; iter++)
    fib[iter] = fib[iter-1] + fib[iter - 2];
   return fib[n];
}

int main(int argc, char **argv) {
  int n;
  cout << "What order Legendre polynomial? ";
  cin >> n;
  
  // if(argc = 1)
  //  int n = atoi(argv[1]);
    int nmax = 6; 
     unsigned long   int f = fib(nmax);
     printf("\n  Recursive tree for n = %d  : %lu  calls %llu  \n",nmax, f, r);

     unsigned long  int iterF = fibIterate(n);
    printf("\n Iterative  return: %lu  Number of calls 1 \n", iterF);

   
    unsigned long  int recF = fibClean(n);
    printf("\n  Recursive return: %lu  Number of calls;   %llu\n", recF,r);
   
    return 1;
}
  
