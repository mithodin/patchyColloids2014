/** @file */

#include <math.h>
#include "distance.h"

/**
 * Calculate the distance between two points
 * 
 * @return The distance between two points
 */
double distance(double x1,double z1,double x2,double z2){
	double dx=x1-x2;
	double dz=z1-z2;
	return sqrt(dx*dx+dz*dz);
}
