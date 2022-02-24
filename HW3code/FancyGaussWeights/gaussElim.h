#ifndef GAUSSELIM_H
#define GAUSSELIM_H
#include <iostream>


double* multiply(double** A,double *b,int N,int r,double val){
	double * ret = new double[N+1];
	if(val == 0.0){
		return ret;
	}
	int i=0;
	for(;i<N;i++){
		*(ret+i) = *(*(A+r)+i) * val;	
	}
	*(ret+N) = *(b+r) * val;
	return ret;
}

int subtract(double** A,double *b,int N,int r1,double * temp){
	int i=0;
        for(;i<N;i++){
                *(*(A+r1)+i) -= *(temp+i);
        }
        *(b+r1) -= *(temp+N);
        return 0;
}

int eliminate(double** A, double *b,int N,int c){
	for(int i=0;i<N;i++){
		if(i != c){
			double * temp = multiply(A,b,N,c,*(*(A+i)+c));
			subtract(A,b,N,i,temp);
		}
	}
return 0;
}

int unitify(double** A, double* b,int N,int r){
	double den = *(*(A+r)+r);
	for(int i=0;i<N;i++){
		*(*(A+r)+i) /= den;
	}
	*(b+r) /= den;
	return 0;
}

//Only works for full matrix.
int gaussianElimination(double** A, double* b, int N)
{
	for(int i=0;i< N;i++)
	{
		unitify(A,b,N,i);
		eliminate(A,b,N,i);
	}
	return 0;
}
#endif
