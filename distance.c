#include <math.h>
#include "distance.h"
#include "parameters.h"

double distance(double x1,double z1,double x2,double z2){
	double dx=x1-x2;
	double dz=z1-z2;
	dx -= round(dx / width) * width;
	return sqrt(dx*dx+dz*dz);
}

double realX(double x){
	return x-floor(x/width)*width;
}

double realZ(double z){
	if ( z < 0 ) return 0;
	if ( z > height ) return height;
	else return z;
}

double realDx(double dx){
	dx -= round(dx / width) * width;
	return dx;
}

double realDz(double dz){
	return dz;
}
