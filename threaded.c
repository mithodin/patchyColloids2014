#include <pthreads.h>
#include "parameters.h"
#include "monte_carlo.h"
#include "statistics.h"
#include "initialize.h"

void printPositions(char[40], Colloid *, double);

void *newThread(void *params){
	Config *c = (Config *)params;
	FILE *out = fopen(c->out,"w");
	Colloid *particles = (Colloid *)malloc(sizeof(Colloid)*(c->N));
	Stats *stat = initStats(c->height);
	initDmax(particles,c,out);
	
	double pacc=monteCarloSteps(particles,c->steps,c,stat);
	fprintf(out,"Overall acceptance rate: %f\n",pacc);

	c->Utot = totalEnergy(particles, &(c->Uext), &(c->Uint));
	fprintf(out,"Total Energy: %f\n",c->Utot);

	printPositions(c->posOut, particles, pacc);
	printStats(c->statOut, particles);

	if(collisions(particles)) fprintf(out,"Collision detected in final state!\n");
	
	*done = 1;

	free(particles);
	free(stat);
	fclose(out);
	free(c);
	pthread_exit(0);
}

void printPositions(char fn[40], Colloid *particles, double paccept){
	FILE *posFile = fopen(fn,"w");
	int i=0;
	fprintf(posFile,"#N = %d\tN1 = %d\tN2 = %d\tU0 = %f\tM1 = %f\tM2 = %f\theight = %f\twidth = %f\tg = %f\tT = %f\tSteps = %d\tAcceptance Rate = %f\n",N,N1,N2,U0,M1,M2,height,width,g,T,steps,paccept);
	fprintf(posFile,"#x\tz\tangle\tspecies (3patch: %d)\n",THREEPATCH);
	for(i=0;i<N;i++){
		fprintf(posFile,"%f\t%f\t%f\t%d\n",particles[i].x,particles[i].z,particles[i].a,particles[i].sp);
	}
	fclose(posFile);
}
