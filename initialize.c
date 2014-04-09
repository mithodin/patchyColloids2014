#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "config.h"
#include "random.h"
#include "parameters.h"
#include "colloid.h"
#include "distance.h"
#include "initialize.h"
#include "statistics.h"
#include "monte_carlo.h"

extern FILE *initFile;

void initRandomly(Colloid *, Config*);
void initBoxed(Colloid *, Config*);
void initFromFile(Colloid *, Config *);

void initParticles(Colloid *particles, Config *c){
	if(c->loadInit){
		initFromFile(particles, c);
	}else if(c->boxed == 0){
		initRandomly(particles, c);
	}else{
		initBoxed(particles,c);
	}
	c->Utot = totalEnergy(particles, c);
	printf("initialization done\n");
}

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
		//Discard first two lines.
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
	}
}

void initRandomly(Colloid *particles, Config *c){
	printf("Initalizing randomly\n");
	int i=0;
	Colloid *this = NULL;
	while(i<c->N){
		pthread_mutex_lock( &mtxRandom );
		particles[i].x = dsfmt_genrand_close_open(&randState)*c->width;
		particles[i].z = dsfmt_genrand_close_open(&randState)*(c->height-1)+0.5; //For not bumping into walls
		particles[i].a = dsfmt_genrand_close_open(&randState)*2*M_PI;
		pthread_mutex_unlock( &mtxRandom );
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
	printf("Initalizing boxed\n");
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
		pthread_mutex_lock( &mtxRandom );
		particles[i].x = dsfmt_genrand_close_open(&randState)*c->width;
		particles[i].a = dsfmt_genrand_close_open(&randState)*2.0*M_PI;
		particles[i].z = dsfmt_genrand_close_open(&randState)*(c->height*box-1)+0.5+offset;
		pthread_mutex_unlock( &mtxRandom );
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
