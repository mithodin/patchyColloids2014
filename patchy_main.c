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
const double delta = 0.119; //Patch diameter
const double sigma = 1.0; //Colloid diameter

double Utot = 0; //Total Energy

void printPositions(void);

int main(void){ 	//This is only for testing so far.
	parameters=getParams();
	if(parameters==NULL) return 1;
	
	loadParams(parameters);
	printf("N = %d\nN1 = %d\nN2 = %d\nU0 = %f\nM1 = %f\nM2 = %f\nheight = %f\nwidth = %f\nT = %f\n",N,N1,N2,U0,M1,M2,height,width,T);

	long int seed=random_seed();
	init_genrand((unsigned long)seed);

	particles=(Colloid *)malloc(sizeof(Colloid)*N);
	initParticles(particles);
	printf("Particles initialized\n");
	//printColloidsSortedX(particles);
	//printColloidsSortedZ(particles);
	Utot = totalEnergy(particles);
	printf("Total Energy: %f\n",Utot);

	initDmax(particles);

	double pacc=monteCarloSteps(particles,50000);
	printf("Overall acceptance rate: %f\n",pacc);
	
	Utot = totalEnergy(particles);
	printf("Total Energy: %f\n",Utot);
	printPositions();

	if(collisions(particles)) printf("Collision detected!\n");
	return 0;
}

void printPositions(void){
	FILE *posFile = fopen("positions.dat","w");
	int i=0;
	fprintf(posFile,"#x\tz\tangle\tspecies (3patch: %d)\n",THREEPATCH);
	for(i=0;i<N;i++){
		fprintf(posFile,"%f\t%f\t%f\t%d\n",particles[i].x,particles[i].z,particles[i].a,particles[i].sp);
	}
	fclose(posFile);
}
