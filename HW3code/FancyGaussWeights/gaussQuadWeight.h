#include <iostream>
#include <vector>
#include <cmath>
#include "gaussElim.h"

#define TOL 1.0E-15
using namespace std;



// evaluate a polynomial of order 'n'
double evalPolynomial(double x, double* a, int n)
{
  if (n == 0) { return a[0]; }
  double sum = a[0];
  for(int i=1;i<n+1;i++){
    sum += *(a+i) * pow(x,i);
  }
  return sum;
}
// evaluate the derivative of a polynomial of order 'n'
double evalPolyDeriv(double x, double* a, int n)
{
  if (n == 0) { return 0; }
  if (n == 1) { return a[1]; }
  double sum = a[1];
  for(int i=2;i<n+1;i++){
    sum += *(a+i) * i * pow(x,i-1);
  }
  return sum;
}

double getRootsUsingNewton(double * a,double guess,int n){
	double value = guess;
	int count = 0;
	do{
		count++;
		value = value - (evalPolynomial(value,a,n)/evalPolyDeriv(value,a,n));
	}while(value < TOL && count < 100000);
	return value;
}

vector<double> mul(vector<double> v1, vector<double> v2){
        vector<double> ret(v1.size() + v2.size() -1,0);
        for(int i=0;i<v1.size();i++){
                for(int j=0;j<v2.size();j++){
                        ret[i+j] = ret.at(i+j) + v1.at(i)*v2.at(j);
                }
        }
        return ret;
}

vector<double> sub(vector<double> v1, vector<double> v2){
        int min_size = v1.size()<v2.size() ? v1.size() : v2.size();
        int max_size = v1.size()>v2.size() ? v1.size() : v2.size();
        vector<double> ret(max_size,0);
        int i =0;
        for(;i<min_size ;i++){
                ret[i] =  v1.at(i) - v2.at(i);
        }
        int to_mul  = 1;
        vector<double> rem;
        if(v2.size()>v1.size()){
                rem = v2;
                to_mul = -1;
        } else {
                rem = v1;
        }
        for(;i<rem.size();i++){
                ret[i] = to_mul*rem.at(i);
	}
        return ret;
}

void vecToArr(vector<double> v, double * a){
        for(int i =0;i<v.size();i++){
                *(a+i) = v[i];
        }
}

vector<double> legendre(double n){
        vector<double> ret;
        if(n == 0){
                ret.push_back(1);
        }
        else if(n == 1){
                ret.push_back(0);
                ret.push_back(1);
	}else {
                vector<double> v1;
                v1.push_back(0);
                double temp = ((2*n)-1)/n;
                v1.push_back(((2*n)-1)/n);
                vector<double> v2 = legendre(n-1);
                vector<double> v_res = mul(v1,v2);
                vector<double> v3;
                v3.push_back((n-1)/n);
                vector<double> v4 = legendre(n-2);
                vector<double> v1_res = mul(v3,v4);
                ret = sub(v_res,v1_res);
        }
        return ret;
}

int getLegendreCoeff(double * a, double n){
        if(n<0){
                cout << "Only inetrger values supported" << endl;
                return 0;
        }
        if(a == NULL){
                cout << "the array is NULL" << endl;
                return 0;
        }
        vector<double> val = legendre(n);
        vecToArr(val,a);
        return 1;
}

int getLegendreZero(double* zero, double* a, double n)
{
  if (zero == 0) { return -1; } // error out if not allocated
  if (a == 0) { return -1; } // error out if not allocated
  for(int i=n,j=0;i>0;i--,j++){
    double g1 = cos((((4*i)-1)/((4*n)+2)) * 3.14159265);
    double g2 = 1 - (1/(8*pow(n,2))) + (1/(8*pow(n,3)));
    double guess = g1 * g2;
    *(zero+j) = getRootsUsingNewton(a,guess,n);
  }
  return 0;
}

int getGaussQuadXW(double * x_i,double * w_i,int order){

  double * leg_coeffs = new double[order+1];
  getLegendreCoeff(leg_coeffs,order);
  getLegendreZero(x_i, leg_coeffs, order); 
  
  double **A = new double*[order];
  for(int i=0;i<order;i++)
    *(A+i) = new double[order];
  
    
  for(int i=0;i<order;i++)
    for(int j=0;j<order;j++)
      *(*(A+i)+j) = pow(x_i[j],i);
  
  for(int i=1,j=0;i<order+1;i++,j++)
    *(w_i+j) = (1 + pow(-1,i-1))/i;

  gaussianElimination(A,w_i,order);
  return 0;
}
