#include <math.h>
#include "distance.h"
#include "parameters.h"

double distance(double x1,double z1,double x2,double z2){
	double dx=x1-x2;
	double dz=z1-z2;
	dx -= round(dx / width) * width;
	dz -= round(dz / height) * height;
	return sqrt(dx*dx+dz*dz);
}

double realX(double x){
	return x-floor(x/width)*width;
}

double realZ(double z){
	return z-floor(z/height)*height;
}
