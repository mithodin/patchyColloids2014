#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mt19937ar.h"
#include "parameters.h"
#include "colloid.h"
#include "distance.h"
#include "initialize.h"
#include "monte_carlo.h"

extern FILE *initFile;

void initParticles(Colloid *particles){
	if(initFile){
		initFromFile(particles);
	}else{
		initRandomly(particles);
	}
	makePeriodicX(particles);
	Utot = totalEnergy(particles, &Uext, &Uint);
}

void initRandomly(Colloid *particles){
	int i=0;
	Colloid *this = NULL;
	while(i<N){
		particles[i].x = genrand_real2()*width;
		particles[i].z = genrand_real2()*(height-1)+0.5; //For not bumping into walls
		particles[i].a = genrand_real2()*2*M_PI;
		if(noCollision(i,particles)){
			this = &particles[i];
			if(i<N1){
				newColloid(THREEPATCH,this);
			}else{
				newColloid(TWOPATCH,this);
			}
			if(i>0){
				insertSortedX(&particles[i-1], this);
				insertSortedZ(&particles[i-1], this);
			}
			i+=1;
		}
	}
}

void initFromFile(Colloid *particles){
	//not doing this!
}

int noCollision(int i, Colloid *particles){
	double x = particles[i].x;
	double z = particles[i].z;
	double xj = 0, zj = 0;
	int j=0;
	extPotential(&particles[i],&j);
	if( j ) return 0;
	for(j=0;j<i;j++){
		xj = particles[j].x;
		zj = particles[j].z;
		if(distance(x,z,xj,zj) < 1.0) return 0;
	}
	return 1;
}
