#include <math.h>
#include <stdio.h>
#include "colloid.h"
#include "monte_carlo.h"
#include "parameters.h"
#include "distance.h"

double dmax = 1.0; //completely random numbers
double amax = 2.0/3.0*M_PI;

double extPotential(Colloid *c){
	return 0;
}

double pairPotential(Colloid *particle){
	Colloid *partner = particle;
	double x = (*particle).x;
	double z = (*particle).z;
	double u = 0;
	while( fabs(realDx(x - (*(*partner).left).x)) <= sigma+delta ){
		partner = (*partner).left;
		if( fabs(realDz(z - (*partner).z)) <= sigma+delta ){
			u += pairInteraction(particle,partner);
		}
	}
	partner = particle;
	while( fabs(realDx(x - (*(*partner).right).x)) <= sigma+delta ){
		partner = (*partner).right;
		if( fabs(realDz(z - (*partner).z)) <= sigma+delta ){
			u += pairInteraction(particle,partner);
		}
	}
	return u;
}

double totalEnergy(Colloid *carray){ //Give an Array here!
	double utot = 0;
	double uext = 0;
	double uint = 0;
	int i = 0;
	for(i = 0;i < N; ++i){
		uext += extPotential(&carray[i]);
		uint += pairPotential(&carray[i]);
	}
	utot = uext + uint/2.0;
	return utot;
}

void monteCarloStep(Colloid *carray){
	int i=0;
	for(i = 0; i < N; ++i){
		//TODO
	}
}

void monteCarloSteps(Colloid *carray, int howmany){
	for(; howmany > 0; --howmany){
		monteCarloStep(carray);
	}
}

int accept(double du){
	if( du < 0 ) return 1;
	else{
		double paccept = exp(-du/T);
		//TODO
	}
}
