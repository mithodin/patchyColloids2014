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
#include "headers.h" 	//Global headers file
#include "colloid.h"	//Define what is a colloid?

config_t *parameters;	//holds the parameters for our simulation
Colloid *particles;	//holds all the particles

/* No need to use these so far:
Colloid *pSortedX;	//particles sorted by x coordinate
Colloid *pSortedZ;	//particles sorted by z coordinate
*/

int N;
int N1;
int N2;
double U0;
const double M1 = 1;
double M2;

int main(void){ 	//This is only for testing so far.
	parameters=getParams();
	if(parameters==NULL) return 1;
	
	loadParams(parameters);
	printf("N = %d\nN1 = %d\nN2 = %d\nU0 = %f\nM1 = %f\nM2 = %f\n",N,N1,N2,U0,M1,M2);
	return 0;
}
