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

config_t *parameters;	//holds the parameters for our simulation
Colloid *particles;	//holds all the particles

int N;
int N1;
int N2;
double U0;
const double M1 = 1;
double M2;
double height;
double width;

int main(void){ 	//This is only for testing so far.
	parameters=getParams();
	if(parameters==NULL) return 1;
	
	loadParams(parameters);
	printf("N = %d\nN1 = %d\nN2 = %d\nU0 = %f\nM1 = %f\nM2 = %f\nheight = %f\nwidth = %f\n",N,N1,N2,U0,M1,M2,height,width);

	long int seed=random_seed();
	init_genrand((unsigned long)seed);

	particles=(Colloid *)malloc(sizeof(Colloid)*N);
	initParticles(particles);
	printf("Particles initialized\n");
	int i=0;
	while(i<N){
		printf("%f ",particles[i++].x);
		fflush(stdout);
	}
	printf("\n");
	printColloidsSortedX(particles);
	return 0;
}
