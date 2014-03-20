/*************************************************
 * Main Simulation file for patchy colloids
 * in a gravitational field
 *
 * Bachelor's Thesis by Lucas Treffenstädt
 *************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <libconfig.h> 	//We are using this to load our parameters from a file.

#include "mt19937ar.h" 	//Generating random numbers
#include "colloid.h" 	//What is a colloid?
#include "initialize.h"	//Initialization
#include "load_config.h" //Loading the configuration
#include "monte_carlo.h"
#include "statistics.h"

config_t *parameters;	//holds the parameters for our simulation
Colloid *particles;	//holds all the particles

int N;
int N1;
int N2;
const double M1 = 1;
double M2;
double height;
double width;
double T;
const double U0 = 1;
const double delta = 0.11965683746373795115; //Patch diameter
//const double delta = 0.2;
const double sigma = 1.0; //Colloid diameter

FILE *initFile = NULL;

double Utot = 0; //Total Energy
double Uint = 0;
double Uext = 0;

extern double *rho1,*rho2,*f1,*f2;
extern double dmax,amax,simRate;
extern bool ineq;

int steps;

char fn[40];

void printPositions(double);
void reset(void);

int main(void){ 	//This is only for testing so far.
	parameters=getParams();
	if(parameters==NULL){
		printf("Error loading parameters.cfg\n");
		return 1;
	}
	
	long int seed=random_seed();
	init_genrand((unsigned long)seed);

	int params = 0;
	while( (params = loadParams()) ){;
		printf("Round %d\n",params);
		printf("N = %d\nN1 = %d\nN2 = %d\nU0 = %f\nM1 = %f\nM2 = %f\nheight = %f\nwidth = %f\nT = %f\nSteps = %d\n",N,N1,N2,U0,M1,M2,height,width,T,steps);

		particles=(Colloid *)malloc(sizeof(Colloid)*N);
		initParticles(particles);
		printf("Particles initialized\n");
		printf("Total Energy: %f\n",Utot);

		initStats(height);

		initDmax(particles);

		double pacc=monteCarloSteps(particles,steps);
		printf("Overall acceptance rate: %f\n",pacc);

		Utot = totalEnergy(particles, &Uext, &Uint);
		printf("Total Energy: %f\n",Utot);
		printPositions(pacc);

		printStats(particles);

		if(collisions(particles)) printf("Collision detected!\n");
		reset();
	}
	printf("So long and thanks for all the fish...\n");
	return 0;
}

void printPositions(double paccept){
	FILE *posFile = fopen(fn,"w");
	int i=0;
	fprintf(posFile,"#N = %d\tN1 = %d\tN2 = %d\tU0 = %f\tM1 = %f\tM2 = %f\theight = %f\twidth = %f\tT = %f\tSteps = %d\tAcceptance Rate = %f\n",N,N1,N2,U0,M1,M2,height,width,T,steps,paccept);
	fprintf(posFile,"#x\tz\tangle\tspecies (3patch: %d)\n",THREEPATCH);
	for(i=0;i<N;i++){
		fprintf(posFile,"%f\t%f\t%f\t%d\n",particles[i].x,particles[i].z,particles[i].a,particles[i].sp);
	}
	fclose(posFile);
}

void reset(void){
	free(particles);
	free(rho1);
	free(rho2);
	free(f1);
	free(f2);
	dmax = 1.0;
	amax = 0.1*2.0/3.0*M_PI;
	simRate = 0;
	ineq = false;
}
