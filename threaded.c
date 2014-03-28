#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include "config.h"
#include "colloid.h"
#include "statistics.h"
#include "parameters.h"
#include "monte_carlo.h"
#include "initialize.h"

void printPositions(Colloid *, double, Config *);

void *newThread(void *params){
	printf("I'm alive!\n");
	Config *c = (Config *)params;
	FILE *out = fopen(c->out,"w");
	if( !out ){
		printf("Error opening output file\n");
		fclose(out);
		pthread_exit((void *)1);
	}else{
		fprintf(out,"Thread alive\n");
		fflush(out);
	}
	Colloid *particles = (Colloid *)malloc(sizeof(Colloid)*(c->N));
	fprintf(out,"Initializing particles\n");
	fflush(out);
	initParticles(particles,c);
	fprintf(out,"Particles initialized!\n");
	fflush(out);
	fclose(out);

	Stats *stat = initStats(c->height/2);
	initDmax(particles,c);
	
	out = fopen(c->out,"a");
	fprintf(out,"Running simulation\n");
	fflush(out);
	double pacc=monteCarloSteps(particles,c->steps,c,stat,out);
	fprintf(out,"Overall acceptance rate: %f\n",pacc);
	fflush(out);

	c->Utot = totalEnergy(particles, c);
	fprintf(out,"Total Energy: %f\n",c->Utot);
	fflush(out);

	printPositions(particles, pacc,c);
	printStats(particles, c->height, stat, c->statOut);

	if(collisions(particles,c)) fprintf(out,"Collision detected in final state!\n");
	
	fprintf(out,"I'm done!\n");
	*(c->done) = 1;

	free(particles);
	free(stat);
	fclose(out);
	free(c);
	pthread_exit((void *)0);
}

void printPositions(Colloid *particles, double paccept, Config *c){
	FILE *posFile = fopen(c->posOut,"w");
	int i=0;
	fprintf(posFile,"#N = %d\tN1 = %d\tN2 = %d\tU0 = %f\tM1 = %f\tM2 = %f\theight = %f\twidth = %f\tg = %f\tT = %f\tSteps = %d\tAcceptance Rate = %f\n",c->N,c->N1,c->N2,U0,M1,c->M2,c->height,c->width,c->g,c->T,c->steps,paccept);
	fprintf(posFile,"#x\tz\tangle\tspecies (3patch: %d)\n",THREEPATCH);
	for(i=0;i<c->N;i++){
		fprintf(posFile,"%f\t%f\t%f\t%d\n",particles[i].x,particles[i].z,particles[i].a,particles[i].sp);
	}
	fclose(posFile);
}
