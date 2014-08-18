/** @file */

#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include "config.h"
#include "colloid.h"
#include "statistics.h"
#include "parameters.h"
#include "monte_carlo.h"
#include "initialize.h"
#import "dSFMT/dSFMT.h"

void printPositions(Colloid *, double, Config *);
int getRandomSeed(void);

bool initdmax = true; /**< initialize amax and dmax? */

/**
 * Spawns a new thread with given parameters
 *
 * @param params Pointer to Config struct. Make sure to initialize correctly (see load_config.c)
 */
void *newThread(void *params){
	//printf("I'm alive!\n");
	Config *c = (Config *)params;
	dsfmt_init_gen_rand(&(c->myrand),getRandomSeed()); //initialize the rng
	FILE *out = fopen(c->out,"w");
	if( !out ){
		printf("Error opening output file\n");
		fclose(out);
		pthread_exit((void *)1);
	}else{
		fprintf(out,"Thread alive\n");
		fprintf(out,"#N = %d\nN1 = %d\nN2 = %d\nU0 = %f\nM1 = %f\nM2 = %f\nheight = %f\nwidth = %f\ng = %f\tT = %f\nSteps = %d\n",c->N,c->N1,c->N2,U0,M1,c->M2,c->height,c->width,c->g,c->T,c->steps);
		fflush(out);
	}
	Colloid *particles = (Colloid *)malloc(sizeof(Colloid)*(c->N));
	fprintf(out,"Initializing particles\n");
	fflush(out);
	initParticles(particles,c); //initialize the particle positions as defined by config
	fprintf(out,"Particles initialized!\n");
	fflush(out);
	fclose(out);
	Stats *stat = initStats((int)(c->height/2.0));

	//if we want to skip the initialization part for debugging purposes
	#ifdef PM_DEBUG_NOINIT
	c->dmax = 5e-2;
	c->amax = 5e-1;
	#else
	if(initdmax){
		initDmax(particles,c); //initialize maximum displacement and rotation
	}
	#endif
	
	out = fopen(c->out,"a");
	fprintf(out,"Running simulation\n");
	fflush(out);
	double pacc=monteCarloSteps(particles,c->steps,c,stat,out); //run the main simulation
	fprintf(out,"Overall acceptance rate: %f\n",pacc);
	fflush(out);

	c->Utot = totalEnergy(particles, c); //calculate the total energy in the final state
	fprintf(out,"Total Energy: %f\n",c->Utot);
	fflush(out);

	printPositions(particles, pacc,c); //print the final configuration to file
	printStats(stat, c); //print the statistics (density,bond saturation) to file

	if(collisions(particles,c)) fprintf(out,"Collision detected in final state!\n");
	
	fprintf(out,"I'm done!\n");
	*(c->done) = 1; //signal the main process that I'm ready to join

	//Free all our resources.
	free(particles);
	free(stat->rho1);
	free(stat->rho2);
	free(stat->f1);
	free(stat->f2);
	free(stat);
	fclose(out);
	free(c);
	pthread_exit((void *)0);
}

/**
 * Print positions of a given state to file
 *
 * @param particles Pointer to array of Colloids of length c->N
 * @param paccept Acceptance rate to include in the header
 * @param c Pointer to Config struct
 */
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

/**
 * Get a random int to use as seed
 *
 * If /dev/random is not available or fails to read, uses timeofday to get a "random" number
 *
 * @return Returns a more or less random int
 */
int getRandomSeed(void){
	FILE *rnd = fopen("/dev/random","r");
	int seed;
	if( rnd ){
		if(fread(&seed,sizeof(int),1,rnd) != 1){
			struct timeval t;
			gettimeofday(&t,NULL);
			seed=((int)t.tv_usec)^((int)t.tv_sec);
		}
	}else{
		struct timeval t;
		gettimeofday(&t,NULL);
		seed=((int)t.tv_usec)^((int)t.tv_sec);
	}
	return seed;
}
