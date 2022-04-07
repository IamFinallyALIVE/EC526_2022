/**
 * Author: Michael Graziano mjgrazia@bu.edu
 * Author: Sahan Bandara sahanb@bu.edu
 */

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <chrono>
#include <vector>
#include <cmath>

#define R_cut 2.0
#define R_skin 4.0

using namespace std;

// Function Prototypes:
int main(int argc, char **argv );
void setPositions (int np, int nd, double L, double *pos);
double sign (double a, double b); 
void createNeighborhood (int np, int nd, double *pos, vector<double> *neighbor);
void calcAccel (int np, int nd, double L, double mass, double *pos, double *acc, double &Epot);
void boundCond (int np, int nd, double L, double *pos);
void calcVerlet (int np, int nd, double mass, double dt, double L, double &Ekin, double &Epot,
  double *pos, double *vel, double *acc);

int main(int argc, char **argv) {
    /**
     * Main function of the Molecular Dynamics Simulator
     */

    // Variable Declartions
    double *pos;              // Holds the position of the particles
    vector<double> *neighbor; // Holds the neighborhood table for the simulation
    double *vel;              // Holds the velocity of the particles
    double *acc;              // Holds the accleration of the particles
    double L = 4096;           // Size of the simulation area
    double Ekin, Epot;        // Holds the kinetic & potential energy of the system.
    double r_ij;              // Holds the distance between pairs of particles
    double mass = 1.0;        // Holds the mass of each atom
    double dt = 1E-12;        // Holds the time step
    double steps = 5000;      // Holds the maximum number of time steps to do
    int np = 2048;               // Holds the number of particles
    int nd = 2;               // Holds the number of dimensions
    int i, j, k, dim;         // General iterators
    
    // Chrono Time Variables
    chrono::time_point<chrono::steady_clock> begin_time, end_time;
    chrono::duration<double> difference_in_time;
    double difference_in_seconds;

    // Initialize memory
    pos = new double[nd*np];
    vel = new double[nd*np];
    acc = new double[nd*np];
    neighbor = new vector<double>[np];

    // Initialize Ekin
    Ekin = 0;

    // Initialize velocities
    for (i = 0; i < np; i++) {
        for (dim = 0; dim < nd; dim++) {
            vel[dim+i*nd] = 0.0;
        }
    }
    
    // Initialize positions
    setPositions(np, nd, L, pos);

//    // Create neighbor table
//    createNeighborhood(np, nd, pos, neighbor);

//    // Checking neighbor table creation
//    cout << "Neighborhood Table: " << endl;
//    for (i = 0; i < np; i++) {
//        cout << "Particle #" << i << ": ";
//        for (j = 0; j < neighbor[i].size(); j++) {
//            cout << neighbor[i][j];
//            if (j == neighbor[i].size()-1) { cout << endl; }
//            else { cout << ", "; }
//        }
//    }
            

    // Initialize potential energy
    calcAccel(np, nd, L, mass, pos, acc, Epot);

    // Initialize acceleration to zero
    for (i = 0; i < np; i++) {
        for (dim = 0; dim < nd; dim++) {
            acc[dim+i*nd] = 0.0;
        }
    }
    
    begin_time = chrono::steady_clock::now();

//    cout << scientific; //<< setprecision(3);
    // Perform Verlet Velocity Calculations
    for (i = 0; i < steps; i++) {
//        if ((i+1)%100 == 0 || i == 0 || i == steps-1) {
//            cout << "Step #" << i+1 << " --------------------" << endl;
//            cout << "Potential Energy: " << Epot << "\tKinetic Energy: " << Ekin << "\tTotal Energy: " << Epot + Ekin << endl;
//            cout << "Particle | Position (X) | Position (Y) | Velocity (X) | Velocity (Y) | Acceler. (X) | Acceler. (Y) |" << endl;
//            cout << "Particle | Position (X) | Velocity (X) | Acceler. (X) |" << endl;
//            for (j = 0; j < np; j++) {
//                cout << resetiosflags(ios::scientific) << setw(8) << j+1 << " |";
//                cout << scientific;
//                for (k = 0; k < 3; k++) {
//                    for (dim = 0; dim < nd; dim++) {
//                        if (k == 0) { // Position
//                            cout << setw(13) << pos[dim+j*nd];
//                        } else if (k == 1) {
//                            cout << setw(13) << vel[dim+j*nd];
//                        } else {
//                            cout << setw(13) << acc[dim+j*nd];
//                        }
//                        cout << " |";
//                    }
//                }
//                cout << endl;
//            }
//        }
        calcVerlet (np, nd, mass, dt, L, Ekin, Epot, pos, vel, acc);
    }

    end_time = chrono::steady_clock::now();
    difference_in_time = end_time - begin_time;
    difference_in_seconds = difference_in_time.count();

    cout << np << "\t" << L << "\t" 
         << scientific << setprecision(3) << dt << "\t" 
         << resetiosflags(ios::scientific) << setprecision(15) << difference_in_seconds << endl;

//    // For the time being assuming 1 Dimension. Will need to update for multiple dimensions
//    double avg_v = 0;
//    for (i = 0; i < np; i++) {
//        for (dim = 0; dim < nd; dim++) {
//            avg_v += abs(vel[dim+i*nd]);
//        }
//    }
//    avg_v /= np;
//    cout << "Average Velocity = " << avg_v << endl;
//    cout << "Neighbor Update Time can be < " << (R_skin-R_cut)/(avg_v*dt) << endl;
//            
//    cout << "TEST COMPLETED!" << endl;


    // Clean memory
    delete [] pos;
    delete [] vel;
    delete [] acc;

    return 0;
}

void setPositions (int np, int nd, double L, double *pos) {
    /**
     * Used to randomly intialize the positions of the particles in the system.
     * Code adapted from csm_ch08.pdf.
     */

    // Variable Declarations
    double min_r = pow(2.0, 1.0/3.0); // Used to prevent overlap
    double dist_diff[nd];             // Used to hold the difference in distance
    double overlap_check;             // Used to hold the computation for overlap checking
    bool overlap;                     // Used to indicate overlap occured
    int i, j, dim;                    // General iterators
    
    for (i = 0; i < np; i++) { // Scans through all particles
        do {
            overlap = false;
            for (dim = 0; dim < nd; dim++) { // Scans through all dimensions of particle
                pos[dim+i*nd] = L*static_cast<double>(rand())/static_cast<double>(RAND_MAX);
            }
            j = 0;
            while (j < i && !overlap) {
                overlap_check = 0;
                for (dim = 0; dim < nd; dim++) {
                    overlap_check += (pos[dim+i*nd]-pos[dim+j*nd])
                                     * (pos[dim+i*nd]-pos[dim+j*nd]);
                }
                if (overlap_check < min_r) { overlap = true; }
                j++;
            }
        } while (overlap);
    }
}

void createNeighborhood (int np, int nd, double *pos, vector<double> *neighbor) {
    /**
     * Used to create the neighborhood table for the simulation.
     */

    // Variable Declarations
    int i, j, dim;               // General iterators
    double diff_dist;            // Difference in distance
    bool in_range;               // Indicates if a particle is a neighbor
    
    for (i = 0; i < np; i++) {
        for (j = 0; j < np; j++) {
//        for(j = i+1; j < np; j++) {
            if (j == i) { continue; } // Can't pair the particle with itself
            else { // Check distances
                for (dim = 0; dim < nd; dim++) {
                    diff_dist = pos[dim+i*nd] - pos[dim+j*nd];
                    in_range = (abs(diff_dist) <= R_skin);
                    if (!in_range) { break; }
                }
                if (in_range) { neighbor[i].push_back(j); }
            }
        }
    }
}


double sign (double a, double b) {
    /**
     * Returns the value of a with the sign of b. Follows the same functionality as
     * that described in the Fortran sign function.
     * https://gcc.gnu.org/onlinedocs/gfortran/SIGN.html
     */

    // Variable Declaration:
    double output; // Used to hold the output for the function.

    if (b >= 0) { output = abs(a); }
    else { output = -abs(a); }

    return output;
}

void calcAccel (int np, int nd, double L, double mass, double *pos, double *acc, double &Epot) {
    /**
     * Used to calculate the forces (aka acceleration) based on the interations
     * between particles.
     */

    // Variable Declarations
    double dist_diff[nd];        // Holds the difference in distance
    double r2, r2_over, r6_over; // Holds the r^2, r2_over, r6_over for force calculations
    double fp;                   // Holds the coefficient for the force equation  
    double finst[nd];            // Holds the instantanous force (in each direction)
    bool in_range;               // Holds the flag indicating the particle is in range
    int i, j, dim;               // General iterators
    
    // Zero out force & Epot
    Epot = 0.0;
    for (i = 0; i < np; i++) {
        for (dim = 0; dim < nd; dim++) {
            acc[dim+i*nd] = 0.0;
        }
    }

    // Perform force & Epot calculation based on proximity
    for (i = 0; i < np; i++) {
        for (j = i+1; j < np; j++) {
            r2 = 0.0;
            in_range = true;
            for (dim = 0; dim < nd; dim++) {
                dist_diff[dim] = pos[dim+i*nd] - pos[dim+j*nd];
                // Minimum Image Criterion
                if (abs(dist_diff[dim]) > 0.5*L) { dist_diff[dim] -= sign(L, dist_diff[dim]); }
                if (dist_diff[dim] > R_cut) { in_range = false; }
                r2 += dist_diff[dim] * dist_diff[dim];
            }
            if (in_range) { 
                r2_over = 1.0/r2;
                r6_over = r2_over*r2_over*r2_over;
                fp = 48.0*r6_over*(r6_over-0.5)*r2_over;
                for (dim = 0; dim < nd && in_range; dim++) {
                    finst[dim] = fp * dist_diff[dim];
                    acc[dim+i*nd] += finst[dim]/mass;
                    acc[dim+j*nd] -= finst[dim]/mass;
                }
                Epot += 4.0*(r6_over*r6_over-r6_over);
            }
        }
    }
}

void boundCond (int np, int nd, double L, double *pos) {
    /**
     * Used to reposition particles that exit the simulation area.
     */

    // Variable Declaration
    int i, dim; // General Iterators
    
    for (i = 0; i < np; i++) {
        for (dim = 0; dim < nd; dim++) {
            if (pos[dim+i*nd] < 0) { // Particle has exited the lower boundary.
                pos[dim+i*nd] = fmod(pos[dim+i*nd], L)+L;
            } else if (pos[dim+i*nd] > L) { // Particle has exited the upper boundary.
                pos[dim+i*nd] = fmod(pos[dim+i*nd], L);
            } else { // Particle is right on a boundary or in the simulation region.
                continue;
            }
        }
    }
}

void calcVerlet (int np, int nd, double mass, double dt, double L, double &Ekin, double &Epot,
  double *pos, double *vel, double *acc) {
    /**
     * Used to perform the Verlet velocity calculation (based off the leap frog method)
     */

    // Variable Declarations
    double dt2 = 0.5 * dt; // Holds the half time step for the velocity calculation
    int i, dim;            // General iterators
    
    for (i = 0; i < np; i++) {
        for (dim = 0; dim < nd; dim++) { 
            vel[dim+i*nd] += acc[dim+i*nd] * dt2; 
            pos[dim+i*nd] += vel[dim+i*nd] * dt;
        }
    }
    boundCond (np, nd, L, pos);
    calcAccel (np, nd, L, mass, pos, acc, Epot);
    Ekin = 0.0;
    for (i = 0; i < np; i++) {
        for (dim = 0; dim < nd; dim++) { 
            vel[dim+i*nd] += acc[dim+i*nd] * dt2;
            Ekin += vel[dim+i*nd] * vel[dim+i*nd];
        }
    }
    Ekin = 0.5 * Ekin/np;
}
