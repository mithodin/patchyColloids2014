#include <math.h>
#include "distance.h"
#include "parameters.h"

double distance(double x1,double z1,double x2,double z2){
	double dx=x1-x2;
	double dz=z1-z2;
	dx -= nearbyint(dx / width) * width;
	dz -= nearbyint(dz / height) * height;
	return dx*dx+dz*dz;
}
