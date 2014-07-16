/** @file */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "config.h"
#include "parameters.h"
#include "colloid.h"
#include "distance.h"
#include "initialize.h"
#include "statistics.h"
#include "monte_carlo.h"
#import "dSFMT/dSFMT.h"

void initRandomly(Colloid *, Config*);
void initBoxed(Colloid *, Config*);
void initFromFile(Colloid *, Config *);
void thermalize(Colloid *, Config *);

/**
 * Wrapper function to handle initialization
 *
 * @param particles Array of Colloids (memory must have already been assigned)
 * @param c Configuration struct
 */
void initParticles(Colloid *particles, Config *c){
	if(c->loadInit){
		initFromFile(particles, c);
	}else if(c->boxed == 0){
		initRandomly(particles, c);
		thermalize(particles, c);
	}else{
		initBoxed(particles,c);
		thermalize(particles, c);
	}
	c->Utot = totalEnergy(particles, c);
	printf("initialization done\n");
}

/**
 * Thermalize a random/boxed initial state
 *
 * @param particles Array of Colloids
 * @param c Configuration struct
 */
void thermalize(Colloid *particles, Config *c){
	double g = c->g, T = c->T;
	c->g = 0;
	c->T = 0.5;
	c->dmax = 5e-1;
	c->amax = 1.0;
	monteCarloSteps(particles, 10000, c, NULL, NULL);
	c->g = g;
	c->T = T;
}

/**
 * Initialize the colloids by reading positions from a file
 * There are no checks to ensure that no collisions occur.
 * See line 57 for format.
 *
 * @param particles Array of colloids, length c->N
 * @param c Configuration struct
 */
void initFromFile(Colloid *particles, Config *c){
	printf("Initalizing from file\n");
	FILE *in = fopen(c->initIn,"r");
	if( !in ){
		printf("Initialization file not found or unable to load!\n");
		exit(1);
	}else{
		char discard[400];
		fgets(discard,400,in);
		fgets(discard,400,in);
		//Discard first two lines, they contain comments only.
		int i,s;
		for(i = 0;i<c->N;++i){
			if(fscanf(in,"%lf %lf %lf %d",&(particles[i].x),&(particles[i].z),&(particles[i].a),&s) == 4){;
				switch(s){
					case 0:
						newColloid(TWOPATCH,&particles[i]);
						break;
					case 1:
						newColloid(THREEPATCH,&particles[i]);
						break;
					default:
						printf("Error reading positions from file.\n");
						exit(1);
				}
				if(i>0){
					insertSortedZ(&particles[i-1], &particles[i]);
				}
			}else{
				printf("Error reading positions from file!\n");
				exit(1);
			}
		}
		totalEnergy(particles,c); //We need this to have the correct bonds and energies.
	}
}

/**
 * Initialize the colloids by generating random initial positions
 */
void initRandomly(Colloid *particles, Config *c){
	printf("Initalizing randomly\n");
	int i=0;
	Colloid *this = NULL;
	while(i<c->N){
		particles[i].x = dsfmt_genrand_close_open(&(c->myrand))*c->width;
		particles[i].z = dsfmt_genrand_close_open(&(c->myrand))*(c->height-1)+0.5; //For not bumping into walls
		particles[i].a = dsfmt_genrand_close_open(&(c->myrand))*2*M_PI;
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

/**
 * Initialize the colloids by generating random initial positions in two boxes as defined by c->boxed
 * c->boxed < 0 means three-patch-particles on top
 * c->boxed > 0 means two-patch-particles on top
 * |c->boxed| is always the fraction of height for the lower box
 */
void initBoxed(Colloid *particles, Config *c){
	printf("Initalizing boxed\n");
	int i=0;
	Colloid *this=NULL;
	double box = c->boxed > 0 ? c->boxed : 1+c->boxed;
	double offset = c->boxed > 0 ? 0:-1.0*c->boxed*c->height;
	
	species sp = THREEPATCH;
	while(i<c->N){
		if(i==c->N1){
			sp = TWOPATCH;
			offset = c->boxed < 0 ? 0:box*c->height;
			box = 1-box;
		}
		particles[i].x = dsfmt_genrand_close_open(&(c->myrand))*c->width;
		particles[i].a = dsfmt_genrand_close_open(&(c->myrand))*2.0*M_PI;
		particles[i].z = dsfmt_genrand_close_open(&(c->myrand))*(c->height*box-1)+0.5+offset;
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

/**
 * Check if no collisions with other particles exis
 * 
 * @param i Index of current particle
 * @param particles Array of colloids, being read up to but not including particles[i]
 * @param c Configuration struct
 * @return 1 if there are no collisions, 0 otherwise
 */
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
