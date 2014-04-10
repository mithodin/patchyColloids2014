/**
 * @file statistics.h
 * @brief Defines the Stats struct and the functions to handle it
 */

/**
 * Statistics struct
 */
struct stat{
	long double *rho1; /**< Store the density of species 1 (threepatch) */
	long double *rho2; /**< Store the density of species 2 (twopatch) */
	long double *f1; /**< Store the bonds of species 1 */
	long double *f2; /**< Store the bonds of species 2 */
	int bins; /**< Number of histogram bins */
	long samplingCount; /**< Count how many times the values have been sampled */
};

typedef struct stat Stats;

void printStats(Stats*, Config *);
Stats *initStats(int);
void updateDensity(double, species, Config*, Stats*);
void updateF(double, double, species, Config*, Stats*);
