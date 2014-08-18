/*************************************************
 * Main Simulation file for patchy colloids
 * in a gravitational field
 *
 * Bachelor's Thesis by Lucas Treffenstädt
 *************************************************/
/** @file */
/*! \mainpage Patchy Colloids in a Gravitational Field
 *
 * \section intro Bachelor Thesis by Lucas Treffenstädt
 * This is the documentation of the program patchy_main. Starting at 'Files' is maybe the best.
 * \section build How to build
 * - Install libconfig development files (libconfig-dev package for debian)
 * - Install cmake
 * - Create a directory 'build' within this source directory
 * - Descend into 'build' and execute 'cmake ..'
 * - Then execute 'make'\n
 * \n
 * You should now have an executable 'patchy_main' in the 'bin' folder.
 * \section Usage
 * - Edit the file 'bin/parameters.cfg' (see: Parameters)
 * - cd into 'bin'
 * - Run program using './patchy_main x' (Where x is the number of concurrent threads to run)
 * \n
 * \section params Parameters
 * How to use parameters.cfg
 * \subsection General Syntax
 * &lt;name&gt; = &lt;value&gt;\n
 * where &lt;value&gt; can be a single value (integer, float, string)
 * - float always requires a decimal point! Example: 1.0
 * - integer ist just any integer number. Example: 500
 * - string must be put in single or double quotes. Example: "file.dat"
 * some values can also be arrays. An array starts with a square bracket [ and ends with one ]. Values are comma separated. Example: [1.0, 2.0, 3.0]. Single values are not allowed!
 * \subsection vars Variables
 * - \b N: (integer, array) Total number of particles
 * - \b N1, \b N2: (integer) Number of species 1/2
 * - \b x: (float, array) Composition (molar fraction of species 1)
 * - \b height: (float) Height of the simulation box
 * - \b width: (float) Width of the simulation box
 * - \b T: (float, array) Temperature of the simulation
 * - \b steps: (integer) Number of Monte Carlo steps to run
 * - \b g: (float) Gravitational constant
 * - \b init: (string) Initial configuration
 * 	- "random": random configuration. Will be thermalized
 * 	- "boxed": random configuration, but all species 1 in one part and species 2 in the other part of the box. See "sep". Will be thermalized.
 * 	- "file": load configuration from a file. Will not be thermalized
 * - \b initfile: (string) if init is set to "file", read initial configuration from this file.
 * - \b sep: (float) if init is set to "boxed", this sets the splitting of the simulation box. If positive, species 2 is on top. If negative, species 1 is on top. |sep| is the fraction of height of the lower box.
 * \subsection example Example:
 * <CODE>
 * N = 5000\n
 * x = 0.6\n
 * M2 = -4.12\n
 * height= 200.0\n
 * width= 100.0\n
 * T = 0.1\n
 * steps = 5000\n
 * init = "random"\n
 * g = 0.005\n
 * </CODE>
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <libconfig.h> 	//We are using this to load our parameters from a file.
#include <pthread.h>
#include <unistd.h>

#include "dSFMT/dSFMT.h" 	//Generating random numbers
#include "config.h"
#include "colloid.h" 	//What is a colloid?
#include "statistics.h"
#include "initialize.h"	//Initialization
#include "load_config.h" //Loading the configuration
#include "monte_carlo.h"
#include "threaded.h"
#include "dSFMT/dSFMT.h"

config_t *parameters;	/**< holds the parameters for our simulation */
Colloid *particles;	/**< holds all the particles */

//PARAMS
const double M1 = 1; /**< Mass of species 1 */
const double U0 = 1; /**< Bonding energy */
const double delta = 0.11965683746373795115; /**< Patch diameter */
const double sigma = 1.0; /**< Colloid diameter */
// END PARAMS

int maxThreads = 1; /**< Maximum number of concurrent threads */

int getFreeThread(int,volatile int*, pthread_t*);

/**
 * Main function of the program.
 * @param argc Number of arguments
 * @param argv Array of command line arguments (length argc)
 * @return 0 if program exited successfully, non-zero otherwise
 */
int main(int argc, char** argv){
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
	
	int index = 0;
	int realindex = 0;
	pthread_attr_t at;
	Config **c = (Config **)malloc(sizeof(Config *));
	while( ( index = loadParams(c) ) ){
		realindex=getFreeThread(index,done,threads);
		(*c)->done = &done[realindex];
		done[realindex] = 0;
		sprintf((*c)->out,"%d.stdout",index);

		pthread_attr_init(&at);
		pthread_attr_setdetachstate(&at, PTHREAD_CREATE_JOINABLE);
		int rc = pthread_create(&threads[realindex],&at,newThread,(void *)(*c));
		if( rc ) printf("Error creating thread %d\n",realindex);
		else printf("Thread %d started.\n",realindex);
	}
	int i;
	long result = 0;
	for(i = 0;i<maxThreads;++i){
		pthread_join(threads[i],(void *)result);
		if( result ) printf("Warning: Nonzero exit of thread %d\n",i);
		else printf("thread %d joined and exited normally.\n",i);
	}
	printf("So long and thanks for all the fish...\n");
	return 0;
}

/**
 * Get a new free thread for a simulation
 * @param paramindex Index of the current parameter set
 * @param done Array of ints indicating if a given thread has finished and can be joined
 * @param threads Array of pthreads
 * @return The index of the next free thread
 */
int getFreeThread(int paramindex, volatile int *done, pthread_t *threads){
	if(paramindex > maxThreads){
		while(1){
			int realindex=0;
			long result=0;
			for(realindex = 0;realindex < maxThreads;++realindex){
				if(done[realindex]){
					pthread_join(threads[realindex],(void *)result);
					if( result ){
						printf("Warning: Nonzero exit of thread %d\n",realindex);
					}else{
						printf("Found finished thread\n");
					}
					return realindex;
				}
			}
			sleep(10); //We can afford to wait 10 seconds.
		}
	}else{
		return paramindex-1;
	}
}
