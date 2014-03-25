#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "config.h"
#include "mt19937ar.h"
#include "parameters.h"
#include "colloid.h"
#include "distance.h"
#include "initialize.h"
#include "statistics.h"
#include "monte_carlo.h"

extern FILE *initFile;

void initParticles(Colloid *particles, Config *c){
	initRandomly(particles, c);
	makePeriodicX(particles);
	c->Utot = totalEnergy(particles, c);
}

void initRandomly(Colloid *particles, Config *c){
	int i=0;
	Colloid *this = NULL;
	while(i<c->N){
		particles[i].x = genrand_real2()*c->width;
		particles[i].z = genrand_real2()*(c->height-1)+0.5; //For not bumping into walls
		particles[i].a = genrand_real2()*2*M_PI;
		if(noCollision(i,particles,c)){
			this = &particles[i];
			if(i<c->N1){
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

int noCollision(int i, Colloid *particles, Config *c){
	double x = particles[i].x;
	double z = particles[i].z;
	double xj = 0, zj = 0;
	int j=0;
	extPotential(&particles[i],&j,c);
	if( j ) return 0;
	for(j=0;j<i;j++){
		xj = particles[j].x;
		zj = particles[j].z;
		if(distance(x,z,xj,zj,c) < 1.0) return 0;
	}
	return 1;
}
