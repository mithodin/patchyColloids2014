/**
 * @file config.h
 * @brief This file defines the configuration state for a simulation
 */

#import <stdbool.h>
#import "dSFMT/dSFMT.h"

struct conf{
	volatile int *done; /**< signals the main thread that this subthread is done */
	int N; /**< total number of particles */
	int N1; /**< number of particles in species 1 (three patches) */
	int N2; /**< number of particles in species 2 (two patches) */
	double M2; /**< mass of species 2 */
	double height; /**< height of simulation box */
	double width; /**< width of simulation box */
	double T; /**< temperature in units of \f$\frac{U_0}{k_B}\f$ */
	int steps; /**< number of Monte Carlo steps to run */
	double g; /**< gravitiational constant */
	double dmax; /**< maximum displacement */
	double amax; /**< maximum rotation */
	double boxed; /**< if boxed intialization, give hight of lower box as fraction of height */
	
	double Utot; /**< stores the current total energy */
	double Uint; /**< stores the current internal energy */
	double Uext; /**< stores the current external energy */

	char out[40]; /**< filename for the log messages */
	char statOut[40]; /**< filename for the statistics */
	char posOut[40]; /**< filename for the positions */
	char initIn[50]; /**< filename to read initial positions from */

	bool loadInit; /**< true if inital configuration should be read from file */

	int simRate; /**< stores the number of MC steps per second */
	
	dsfmt_t myrand; /**< stores the RNG state */
};

typedef struct conf Config;
