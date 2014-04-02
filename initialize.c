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

void initRandomly(Colloid *, Config*);
void initBoxed(Colloid *, Config*);

void initParticles(Colloid *particles, Config *c){
	if(c->boxed == 0){
		initRandomly(particles, c);
	}else{
		initBoxed(particles,c);
	}
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
				insertSortedZ(&particles[i-1], this);
			}
			i+=1;
		}
	}
}

void initBoxed(Colloid *particles, Config *c){
	int i=0;
	Colloid *this=NULL;
	double box = c->boxed > 0 ? c->boxed : 1+c->boxed;
	double offset = c->boxed > 0 ? 0:fabs(c->boxed*c->height);
	
	species sp = THREEPATCH;
	while(i<c->N){
		if(i==c->N1){
			sp = TWOPATCH;
			offset = c->boxed < 0 ? 0:box*c->height;
			box = 1-box;
		}
		particles[i].x = genrand_real2()*c->width;
		particles[i].a = genrand_real2()*2.0*M_PI;
		particles[i].z = genrand_real2()*(c->height*box-1 + offset)+0.5;
		if(noCollision(i,particles,c)){
			this = &particles[i];
			newColloid(sp,this);
			if(i>0){
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
		if(distance(x,z,xj,zj) < 1.0) return 0;
	}
	return 1;
}
