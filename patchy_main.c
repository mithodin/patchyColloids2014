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
#include <pthread.h>
#include <unistd.h>

#include "mt19937ar.h" 	//Generating random numbers
#include "config.h"
#include "colloid.h" 	//What is a colloid?
#include "statistics.h"
#include "initialize.h"	//Initialization
#include "load_config.h" //Loading the configuration
#include "monte_carlo.h"
#include "threaded.h"

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
int steps;
// END PARAMS

FILE *initFile = NULL;
char fn[40];
char statFn[40];

int main(int argc, char** argv){
	int maxThreads = 1;
	pthread_t *threads;
	volatile int *done;

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
	Config **c = (Config **)malloc(sizeof(Config *));
	while( ( index = loadParams(c) ) ){;
		printf("N = %d\nN1 = %d\nN2 = %d\nU0 = %f\nM1 = %f\nM2 = %f\nheight = %f\nwidth = %f\ng = %f\nT = %f\nSteps = %d\n",N,N1,N2,U0,M1,M2,height,width,g,T,steps);
		if(index > maxThreads){
			realindex = -1;
			do{
				for(realindex = 0;realindex < maxThreads;++realindex){
					if(done[realindex]){
						pthread_join(threads[realindex],(void *)&result);
						if( result ) printf("Warning: Nonzero exit of thread %d\n",realindex);
						goto freethreadfound;
					}
				}
				sleep(2);
			}while(1);
		}else{
			realindex = index-1;
		}
		freethreadfound:
		(*c)->done = &done[realindex];
		sprintf((*c)->out,"%d.stdout",index);

		pthread_attr_init(&at);
		pthread_attr_setdetachstate(&at, PTHREAD_CREATE_JOINABLE);
		int rc = pthread_create(&threads[realindex],&at,newThread,(void *)c);
		if( rc ) printf("Error creating thread %d\n",realindex);
		else printf("Thread %d started.\n",realindex);
	}
	int i;
	maxThreads = maxThreads > index ? index : maxThreads;
	for(i = 0;i<maxThreads;++i){
		pthread_join(threads[i],(void *)&result);
		if( result ) printf("Warning: Nonzero exit of thread %d\n",i);
	}
	printf("So long and thanks for all the fish...\n");
	return 0;
}
