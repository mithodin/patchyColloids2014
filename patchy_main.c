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
#include <pthreads.h>

#include "mt19937ar.h" 	//Generating random numbers
#include "colloid.h" 	//What is a colloid?
#include "initialize.h"	//Initialization
#include "load_config.h" //Loading the configuration
#include "monte_carlo.h"
#include "statistics.h"
#include "config.h"

config_t *parameters;	//holds the parameters for our simulation
Colloid *particles;	//holds all the particles

//PARAMS
int N;
int N1;
int N2;
const double M1 = 1;
double M2;
double height;
double width;
double T;
const double U0 = 1;
const double delta = 0.11965683746373795115; //Patch diameter
//const double delta = 0.2;
const double sigma = 1.0; //Colloid diameter
double g;
// END PARAMS

FILE *initFile = NULL;
char fn[40];
char statFn[40];

int main(int argc, char* argv){
	int maxThreads = 1;
	pthread_t *threads;
	int *done;

	if ( argc == 2 ){
		maxThreads = atoi(argv[1]);
		if(maxThreads == 0){ ++maxThreads; }
	}else{
		printf("What are you trying to tell me? Running single thread.\n");
	}
	threads = (pthread_t *)calloc(maxThreads,sizeof(pthread_t));
	done = (int *)calloc(maxThreads,sizeof(int));

	parameters=getParams();
	if(parameters==NULL){
		printf("Error loading parameters.cfg\n");
		return 1;
	}
	
	long int seed=random_seed();
	init_genrand((unsigned long)seed);


	int index = 0;
	int realindex = 0;
	int result = 0;
	pthread_attr_t at;
	while( ( index = loadParams() ) ){;
		printf("N = %d\nN1 = %d\nN2 = %d\nU0 = %f\nM1 = %f\nM2 = %f\nheight = %f\nwidth = %f\ng = %f\nT = %f\nSteps = %d\n",N,N1,N2,U0,M1,M2,height,width,g,T,steps);
		Config *c = (Config *)malloc(sizeof(Config));
		c->N = N;
		c->N1 = N1;
		c->N2 = N2;
		c->M2 = M2;
		c->height = height;
		c->width = width;
		c->g = g;
		c->steps = steps;
		c->dmax = 0.5;
		c->amax = 1e-1*2.0/3.0*M_PI;
		sprintf(c->out,"%d.stdout",index);
		
		if(index > maxThreads){
			realindex = -1;
			do{
				for(realindex = 0;realindex < maxThreads;++realindex){
					if(done[realindex]){
						pthread_join(threads[realindex],&result);
						if( result ) printf("Warning: Nonzero exit of thread %d\n",realindex);
						goto freethreadfound;
					}
				}
				nanosleep(2e9); //sleep two seconds
			}while(true);
		}else{
			realindex = index-1;
		}
		freethreadfound:
		c->done = &done[realindex];

		at = pthread_attr_init();
		pthread_attr_setdetachstate(at, PTHREAD_CREATE_JOINABLE);
		int rc = pthread_create(&threads[realindex],at,newThread,(void *)c);
		if( rc ) printf("Error creating thread %d\n",realindex);
	}
	pthread_exit(NULL);
	printf("So long and thanks for all the fish...\n");
	return 0;
}

/*
void reset(void){
	free(particles);
	free(rho1);
	free(rho2);
	free(f1);
	free(f2);
	dmax = 1.0;
	amax = 0.1*2.0/3.0*M_PI;
	simRate = 0;
	ineq = false;
	M1dens = 0;
	M2dens = 0;
}
*/
