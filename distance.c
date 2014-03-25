#include <math.h>
#include "config.h"
#include "distance.h"
#include "parameters.h"

double distance(double x1,double z1,double x2,double z2,Config *c){
	double dx=x1-x2;
	double dz=z1-z2;
	dx -= round(dx / c->width) * c->width;
	return sqrt(dx*dx+dz*dz);
}

double realX(double x, double width){
	return x-floor(x/width)*width;
}

double realZ(double z, double height){
	if ( z < 0 ) return 0;
	if ( z > height ) return height;
	else return z;
}

double realDx(double dx, double width){
	dx -= round(dx / width) * width;
	return dx;
}

double realDz(double dz, double height){
	return dz;
}
