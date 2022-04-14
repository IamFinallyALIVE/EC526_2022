/*  
open boundary condition the particle fly apart.
Need to reduce temperture. Should get a hexagonl crystal

Could put in Hexagon PB with pressure.
*/ 

#include <iostream>
#include <iomanip>
#include <chrono>
#include <vector>
#include <cmath>

using namespace std;

#define PI 3.141592653589793


void setPositions( double *pos,double Radius, int N, int nd);
void calcForce (double *force, double *pos,int N, int nd);
void centerSphere(double *pos, int N, int nd);
bool isOutSide(double *pos,double *reflecvel, double Radius, int N, int nd);


int main(int argc, char **argv )
{
  double *pos;
  double *vel;
  double *force;
  double mass = 1.0;
  int nd = 2;
  int N = 16;
  double R_cut = 2.0;
  double  R_skin = 4.0;

  double radius, angle;
  double dt = 0.05;
  int MaxIter = 1000000;
  int skip =    MaxIter/10;
  
  // Initialize memory
  pos = new double[nd*N];
  vel = new double[nd*N];
  force = new double[nd*N];
  
  //Initialize Position and Velocities

    // Volume per particle is A = Pi Rmax^2/N \simeq 4 sigma^2
  
  double  Radius = 4.0*sqrt(4.0*N/PI);
  cout << "Initial Radius = "<< Radius << endl;
  
  setPositions(pos, Radius, N,  nd);
    
  for(int i = 0; i < N; i++){
    radius = (double)rand()/(double)RAND_MAX;
    angle = 2.0 *PI * (double)rand()/(double)RAND_MAX;
    vel[0 + nd*i] =  0.01*radius*cos(2.0*PI*angle);
    vel[1 + nd*i] =  0.01*radius*sin(2.0*PI*angle);	 
  }
  

#if 0 // check intilization (Should put this printing into separate function call

  calcForce (force, pos,N, nd);
  
  printf("Particle | Position (X) | Position (Y) | Velocity (X) | Velocity (Y) | Force. (X) | force (Y) |\n");    
  for (int j = 0; j < N; j++)
    {
      printf(" %4d   ", j);
      printf("%14.8f       %14.8f       %14.8f       %14.8f     %14.8f       %14.8f\n",
	     pos[0+j*nd],pos[1+j*nd],vel[0+j*nd],vel[1+j*nd],force[0+j*nd],force[1+j*nd]);
      
    }
#endif

  //Interate really don't need start and stop 1/2 steps 

  //Start with half step 

 for (int i = 0; i < N; i++) 
        for (int mu = 0; mu < nd; mu++)  
	  pos[mu+ i*nd] += vel[mu+i*nd] * dt/2.0;       

 for(int iter = 0; iter < MaxIter; iter++)
   {   calcForce (force, pos,N, nd);
     for (int i = 0; i < N; i++) {
       for (int mu = 0; mu < nd; mu++) { //drop  centerSphere
	 //	 if(isOutSide(pos,reflecvel, Radius, N, nd))
	 //  {
	 //   vel[mu+i*nd] = reflecvel[mu+i*nd];
	 //  }
	 //else
	   {
	     vel[mu+i*nd] += force[mu+i*nd] * dt/mass;
	   }
	 pos[mu+i*nd] += vel[mu+i*nd] * dt;
       }
     }
     
#if 0
     if(iter%skip == 0)
       {
	 printf("Particle | Position (X) | Position (Y) | Velocity (X) | Velocity (Y) | Force. (X) | force (Y) |\n");    
	 for (int j = 0; j < N; j++)
	   {
	     printf(" %4d   ", j);
	     printf("%14.8f       %14.8f       %14.8f       %14.8f     %14.8f       %14.8f\n",
		    pos[0+j*nd],pos[1+j*nd],vel[0+j*nd],vel[1+j*nd],force[0+j*nd],force[1+j*nd]);
	     
	   }
    centerSphere(pos, N, nd);
    //Print PE, KE whatever to look at and plot the position of molecules
       }
#endif

    if(iter%skip == 0)  centerSphere(pos, N, nd);
    
     
   }
 
  //Stop

   calcForce (force, pos, N, nd);
     for (int i = 0; i < N; i++) {
       for (int mu = 0; mu < nd; mu++) { 
	 vel[mu+i*nd] += force[mu+i*nd] * dt/mass;
       }
     }
     
 

  
 return 0;
}

void calcForce (double *force, double *pos,int N, int nd)
{
  /*

VLJ(r_ij) = 4 epsilon [ (sigma/r)^(12) - (sigma/r)^6] 
set sigma  = epsioin = 1

F_ij =  -dV(r)/dr  = 48 (1/r)^(13) - 24  (1/r)^(7) 

with r_j = 0  and r=  r_i 

Force on i form j put r_i = 0 so F =  - d/dr \sum_j V[r = r_j] 

  */
  double rij[nd];
  double r2, r2_over, r6_over, fp;
  
  for (int i = 0; i < N; i++)
    for(int mu =0;mu < nd; mu++)
      force[mu + i*nd] =0.0;
  
  for (int i = 0; i < N; i++) {
    // sum over force on i due to j
    for (int j = 1; j < N && i != j; j++) {
      r2 = 0.0;
      for (int mu = 0; mu < nd; mu++) {
	rij[mu] =  pos[mu+i*nd] - pos[mu+j*nd];
	r2 += rij[mu] * rij[mu];
      }
      r2_over = 1.0/r2;
      r6_over = r2_over*r2_over*r2_over;
      fp = 48.0*r6_over*(r6_over-0.5)*r2_over;  // Negative at large distances.
      //   cout << " i = " << i << " j = " << j << endl;
      //  cout << "  r2 = " << r2 << "  r2_over  = " << r2_over << "  r6_over = " << r6_over << "  fp = " << fp << endl;
      for (int mu = 0; mu < nd;  mu++) {
	force[mu + i*nd] += fp * rij[mu];
	//	cout << " mu = " << mu << " rij[mu] = " << rij[mu] << "   force = " << 	force[mu + i*nd] << endl;
      }
    }
    
    //  for (int mu = 0; mu < nd;  mu++) 
    // cout << " mu = " << mu << " i =  " << i <<  "  total force = " <<  force[mu + i*nd] << endl;
  }
}

void centerSphere( double *pos, int N, int nd)
{  double center[nd];
  double radial;
  
  for(int mu = 0; mu < nd; mu++) {
    center[mu] = 0.0;
    for(int i = 0; i < N; i++){
      center[mu] += pos[mu + nd*i];
    }
  }
  
  for(int mu = 0; mu < nd; mu++) {
    center[mu] = center[mu]/(double)N;
    for(int i = 0; i < N; i++) {
      pos[mu + nd*i] -= center[mu];
    }
  }

  // Find Max Radial term;
  double MaxRadial= 0.0;
  double rms = 0.0;
  for(int i = 0; i < N; i++){
    radial = 0.0;
    for(int mu = 0; mu < nd; mu++)   radial +=  pos[mu + nd*i]*pos[mu + nd*i];
    rms +=  radial;
    radial = sqrt(radial);
    if(radial > MaxRadial) MaxRadial = radial;
  }
  
  rms = sqrt(rms/(double) N);
  cout << "  MaxRadial = "<< MaxRadial << "   rms = " << rms  << endl;
}


void setPositions (double *pos, double Radius,  int N, int nd)
{
    /**
     * Used to randomly intialize the positions of the particles in the system.
     * Code adapted from csm_ch08.pdf.
     */
    // Variable Declarations
    double min_r = pow(2.0, 1.0/3.0); // Used to prevent overlap : this is min of potential with sigma = 1
    double overlap_check;             // Used to hold the computation for overlap checking
    bool overlap;                     // Used to indicate overlap occured
    double radius, angle;
   
    
    for (int i = 0; i < N; i++) { // Scans through all particles
        do {
	  overlap = false;

	    radius = (double)rand()/(double)RAND_MAX;
	    angle = 2.0 *PI * (double)rand()/(double)RAND_MAX;
            pos[0 + nd*i] =  radius*Radius*cos(2.0*PI*angle);
            pos[1 + nd*i] =  radius*Radius*sin(2.0*PI*angle);
    
            int j = 0;
            while ( j < i && !overlap) {
                overlap_check = 0;
                for (int mu = 0; mu < nd; mu++) {
                    overlap_check += (pos[mu+i*nd]-pos[mu+j*nd])
                                     * (pos[mu+i*nd]-pos[mu+j*nd]);
                }
                if (overlap_check < min_r) { overlap = true; }
                j++;
            }
        } while (overlap);
    }
}

