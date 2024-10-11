#ifndef H_DECLA
#define H_DECLA

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <sys/stat.h>
#include <sys/types.h>

/*========================= Defining functions =========================*/

//Periodic boundary conditions
#define CLP(A,B)   ( (A) < (0) ? (A+B) : ( (A) > (B) ? (A-B) :A ) )
#define CLPcell(A,B)   ( (A) < (0) ? (A+B) : ( (A) > (B-1) ? (A-B) :A ) )
#define CLPdist(A,B) ((A) <(-0.5*B) ? (A+B) : ((A)>(0.5*B) ? (A-B) : A))
//MIN function
#define MIN(A,B)      ( (A) < (B) ? (A) : (B) )
//SIGN function
#define SIGN(A) ( (A)<(0) ? (-1) : (1))
#define PI 3.14159
#define G 9.8                           //acceleration of gravity (in m.s^-2)
/*========================= Simulation size =========================*/

#define N_PART 130   			        // number of grains
#define L 8					        // length of simulation WARNING : must chose multiple of 2
#define H 50					        // height of simulation
#define NCELLX L                        //number of cell coulmns
#define NCELLY H                        //number of cell rows
#define DENSITE 2.5			            // Must be>2.4 (because 2 + 20% polydispersity)
#define WIDTH grid*R*L		            //width in meters
#define HEIGHT grid*R*H		//height in meters

/*========================= Grain properties =========================*/

#define R 0.001    				        // average particle radius (+-20%)
#define POLYDISP 0.2			        //Dans une population (gros ou petits), on a rmax = (1 + POLYDISP)*R
#define grid 2.5                        //size of initial grid
#define RHO 2500.				        // density of grains
#define MASS (double)(RHO*R*R*R*PI)     			        // mass of a grain of radius R
#define ELL 0.5         //Elastic restitution
#define I (double)(2*MASS*R*R/5)			        // inertia momentum
#define KN (double)(1e3 * MASS*G/R)			            // spring constant, set so that failling of 10d gives d/100 of overlap
#define KT 0.285*KN 			        // tangential spring constant
#define GAMMA -2*log(ELL)*sqrt(MASS*KN)/sqrt(log(ELL)*log(ELL)+PI*PI)				    // dashpot constant
#define GAMMA_PREFACTOR 0.              // tangential gamma
#define MU .25					        // friction VEL_COEFF between particles
#define MU_WALL .25					        // friction VEL_COEFF between particles and walls
#define KR 0	                        // rotational spring constant
#define CR 0	                        // rotational viscosity
#define FORCE_COHESION 0;               // cohesion between grains

//MAGNETIC PARAMETERS
#define MAG_FRACTION 0.                // fraction of magnetic grains
#define NMAG (int)(MAG_FRACTION*(N_PART-L/2))   // number of magnetic grains
#define MAG_MU 0.000025			  // Dipole moment of magnetic grains
#define MAG_NEIGHBOR_CELLS 3			//Number of neighbor cells to check for magnetic forces

/*========================= Physical and simulation constants =========================*/

#define DT (double)(0.05*2*sqrt(MASS/KN)) //Time step
#define SEED 1
#define MU0 4*PI*1e-7 			   // magnetic permeability of free space
#define PREF1 3*MU0/(4*PI)              // prefactor for magnetic force
#define PREF2 MU0/(4*PI)              // prefactor for magnetic torque

/*========================= Output Parameters =========================*/

#define PASTPS_EPS 4000		            // save eps file every PASTPS_EPS
#define PASTPS_EPS_END 2000000			// stop saving eps file after PASTPS_EPS_END
#define NB_ITERATION 2000000           // max number of time steps
#define ENERGY_FREQ 1000                // save energy every ENERGY_FREQ
#define POS_FREQ 100000				   // save positions every POS_FREQ
#define SEDIMENTATION_TIME 50000      // time after which grains are considered to be in sedimentation regime
#define MAGNETIZATION_RAMP_TIME 50000 // time during which the magnetic grains are progressively magnetized
#define SCALE 20000.				        // SCALE size on .eps
#define NMAXcontacts 10                 //Max number of contacts per grain

#define DRAWFORCES 1					//Draw forces on grains
#define MAXWALLS 4					  //Maximum number of walls


#define Y0_PISTON 20                   //initial position of the piston
#define V_PISTON 0.0003                  //velocity of the piston in m/s
/*================================================================================*/

// grains are defined as a structure

typedef  struct grain {

	double x;                           //x coordinate
	double xold;                        //x coordinate at previous timestep
	double dx;                          //x-xold

	double y;                           //y coordinate
	double yold;                        //y coordinate at previous timestep
	double dy;                          //y-xold

	double Oz;                          //Angle of rotation
	double Ozold;                       //Angle at previous timestep
	double dOz;                         //Oz-Ozold

	double fx;                          //horizontal force on slider
	double fy;                          //vertical force on slider
	double Mz;                          //Momentum about z

	//Optional properties for visualization and analysis purposes
	double fnorm;                       //Sum of norms of all forces exerted on grain

	double Ray;                         //Radius of grain
	double mass;                        //mass
	double II;                          //Moment of inertia
	double zp;

	int Nb_Contact;                     //Number of contacts
	int fixed;						 	//fixed grains
	int ismag;
	int highlight;

	int contact[NMAXcontacts];          //List of contacts
	int contactpr[NMAXcontacts];        //List of previous contacts
	int contactwallpr[MAXWALLS];		//List of previous contacts with walls (0 or 1)

	//Contact integrals
	double utijx[NMAXcontacts];
	double utijy[NMAXcontacts];
	double utijx_walls[MAXWALLS];
	double utijy_walls[MAXWALLS];
	double utijxpr[NMAXcontacts];
	double utijypr[NMAXcontacts];
	double utijxpr_walls[MAXWALLS];
	double utijypr_walls[MAXWALLS];
    //double utijOz[NMAXcontacts];
	//double utijOzpr[NMAXcontacts];

} grain;

typedef struct wall {
	double x1, y1; // Starting point of the wall
	double x2, y2; // Ending point of the wall
} wall;

wall walls[MAXWALLS]; // Array to store wall data, assuming a maximum of 4 walls

extern void cell(void);                 //Divide space in cells
extern void createps(char * save_folder, int ieps); //Create .eps file
extern void force_collision(int);       //Compute collision forces
extern void force_ext(void);            //Compute external forces
extern void force_cohesion(void);       //Compute cohesive forces
extern void init_graint(void);          //Initialize grains
extern int init_files(void); //Initialize files
extern void integre_Verlet(void);       //Integrate equations of motion using Verlet's algorithm
extern void integre_VerletStop(void);   //Do not integrate equations of motion
extern void search_collision(void);    //Detect collisions
extern void force_mag(double mag_mu);    		//Compute magnetic forces
extern void force_slider(void);         //Compute forces on slider and update
extern void save_data(void);            //Save data
extern void load_walls(void); // Function to load walls from a file
extern void force_walls(void); // Function to compute forces from walls



int i, j, Nj, somme, bcltps, files_initialized;
int err_count, last_err_count, err;
int tpscoll;
int eps_count;
int HoC[NCELLX+1][NCELLY+1];
int list_link[N_PART+1];
double y_piston, theta, mag_mu, Ec_tot, Ec, fy_piston, fy_left, fy_right;

int p0;

char save_folder[255],energyfilename[255],name_force[255];
FILE *pos_eps,*fp_energy, *fforce;;

grain disk[N_PART+1];

#endif
