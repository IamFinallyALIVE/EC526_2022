// 1D Jacobi MultiGrid

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <string>

#define NN 1024

typedef struct{
    int N;
    int Lmax;
    int size[20];
    double a[20];
    double m;
    double scale[20];
  } param_t;

void relax(double *phi, double *res, int lev, int niter, param_t p);
void proj_res(double *res_c, double *rec_f, double *phi_f, int lev,param_t p);
void inter_add(double *phi_f, double *phi_c, int lev,param_t p);
double GetResRoot(double *phi, double *res, int lev, param_t p);

int main()
{  
  std::string outfile = "MGVALUES_1D.dat";
  FILE *fp = fopen(outfile.c_str(), "w");

  double *phi[20], *res[20];
  param_t p;
  int nlev;
  int i,lev;
  
  //set parameters________________________________________
  p.Lmax = 7; // max number of levels
  p.N = NN;  // MUST BE POWER OF 2
  p.m = 0.01;  
  nlev = 0; // NUMBER OF LEVELS:  nlev = 0 give top level alone 
  if(nlev  > p.Lmax){ 
    printf("ERROR More levels than available in lattice! \n");
    return 0; }
  
  
  printf("\n V cycle for %d by %d lattice with nlev = %d out of max  %d \n", p.N, p.N, nlev, p.Lmax); 
  // initialize arrays__________________________________
  p.size[0] = p.N;
  p.a[0] = 1.0;
  p.scale[0] = 1.0/(2.0 + p.m*p.m);
  
  for(lev = 1;lev< p.Lmax+1; lev++) {
    p.size[lev] = p.size[lev-1]/2;
    p.a[lev] = 2.0 * p.a[lev-1];
    p.scale[lev] = 1.0/(2.0 + p.m*p.m*p.a[lev]*p.a[lev]);
     //p.scale[lev] = 1.0/(2.0 + p.m*p.m);
  }
  
  for(lev = 0;lev< p.Lmax+1; lev++)
    {
      phi[lev] = (double *) malloc(p.size[lev] * sizeof(double));
      res[lev] = (double *) malloc(p.size[lev] * sizeof(double));
      for(i = 0;i< p.size[lev];i++)
	{
	  phi[lev][i] = 0.0;
    res[lev][i] = 0.0;
	};
    }  
  
  //ADD SOURCE
  res[0][p.N/2] = 1.0*p.scale[0];  //unit point source in middle of N by N lattice 
  // iterate to solve_____________________________________
  double resmag = 1.0; // not rescaled.
  int ncycle = 0; 
  int n_per_lev = 10;
  resmag = GetResRoot(phi[0],res[0],0,p);
  //printf("At the %d cycle the mag residue is %g \n",ncycle,resmag);
 
  // while(resmag > 0.00001 && ncycle < 10000)
   while(resmag > 0.000001)
    { 
      ncycle +=1; 
      for(lev = 0;lev<nlev; lev++)   //go down 
	{    
       relax(phi[lev],res[lev],lev, n_per_lev,p); // lev = 1, ..., nlev-1  
        proj_res(res[lev + 1], res[lev], phi[lev], lev,p);    // res[lev+1] += P^dag res[lev]
	}

      for(lev = nlev;lev >= 0; lev--)  //come up
	{ 
  	  relax(phi[lev],res[lev],lev, n_per_lev,p);   // lev = nlev -1, ... 0;
	  if(lev > 0) inter_add(phi[lev-1], phi[lev], lev, p);   // phi[lev-1] += error = P phi[lev] and set phi[lev] = 0;
	}
      resmag = GetResRoot(phi[0],res[0],0,p);
      printf("At the %d cycle the mag residue is %g \n",ncycle,resmag);
    }

  for(int i = 0; i < p.size[0]; i++){
    fprintf(fp,"%.d,    %.10f \n",i,phi[0][i]);
  }

  fclose(fp); // Close the output file
  return 0;
}

void relax(double *phi, double *res, int lev, int niter, param_t p)
{  
  int i,j, x,y;
   int L;
   L  = p.size[lev];  
  //Create Temp T array, set to 0
  double phi_temp[L];

  //Run Jacobi ITerations
  for(i=0; i<niter; i++){

    for(x = 0; x < L; x++){
      // phi_temp[x] = ALPHA * (residual[x] + 0.5 * (phi([x+1] + phi[x-1])) + (1 - ALPHA) * PHI[x]
        phi_temp[x] = 0.5*(res[x] + p.scale[lev] * (phi[(x+1)%L] + phi[(x-1+L)%L])) + 0.5*phi[x]; 
    }

    //Copy It Over
    for(x = 0; x < L; x++){
        phi[x] = phi_temp[x];
    }
  }
  return;    
}

void proj_res(double *res_c, double *res_f, double *phi_f,int lev,param_t p)
{  
  int L, Lc, f_off, c_off, x, y;
  L = p.size[lev];
  double r[L]; // temp residue
  Lc = p.size[lev+1];  // course level
  
  //get residue
  for(x = 0; x < L; x++)
      r[x] = res_f[x] -  phi_f[x]  + p.scale[lev]*(phi_f[(x+1)%L] + phi_f[(x-1+L)%L]);
  
  //project residue
  for(x = 0; x< Lc; x++)
      res_c[x] = 0.5*(r[2*x]  + r[(2*x + 1)%L]);

  return;
}

void inter_add(double *phi_f,double *phi_c,int lev,param_t p)
{  
  int L, Lc, x, y;
  Lc = p.size[lev];  // coarse  level
  L = p.size[lev-1]; 
  
  for(x = 0; x< Lc; x++)
    {
	phi_f[2*x]              += phi_c[x];
	phi_f[(2*x + 1)%L]      += phi_c[x];
    }

      
  //set to zero so phi = error 
  for(x = 0; x< Lc; x++){phi_c[x] = 0.0;}

  return;
}

double GetResRoot(double *phi, double *res, int lev, param_t p)
{ //true residue
  int i, x,y;
  double residue;
  double ResRoot = 0.0;
  int L;
  L  = p.size[lev];
  
  for(x = 0; x < L; x++){
      residue = res[x]/p.scale[lev] - phi[x]/p.scale[lev]  + (phi[(x+1)%L] + phi[(x-1+L)%L]);
      ResRoot += residue*residue; // true residue
    }
  return sqrt(ResRoot);    
}
