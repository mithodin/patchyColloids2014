#include <math.h>
#include <stdio.h>
#include "colloid.h"
#include "monte_carlo.h"
#include "parameters.h"
#include "distance.h"
#include "mt19937ar.h"

double dmax = 1.0; //completely random numbers
double amax = 2.0/3.0*M_PI;

double extPotential(Colloid *c){
	return 0;
}

double pairPotential(Colloid *particle, int *collision){
	Colloid *partner = particle;
	double x = (*particle).x;
	double z = (*particle).z;
	double u = 0;
	while( fabs(realDx(x - (*(*partner).left).x)) <= sigma+delta ){
		partner = (*partner).left;
		if( fabs(realDz(z - (*partner).z)) <= sigma+delta ){
			u += pairInteraction(particle,partner,collision);
		}
	}
	partner = particle;
	while( fabs(realDx(x - (*(*partner).right).x)) <= sigma+delta ){
		partner = (*partner).right;
		if( fabs(realDz(z - (*partner).z)) <= sigma+delta ){
			u += pairInteraction(particle,partner,collision);
		}
	}
	return u;
}

double totalEnergy(Colloid *carray){ //Give an Array here!
	double utot = 0;
	double uext = 0;
	double uint = 0;
	int collision = 0;
	int i = 0;
	for(i = 0;i < N; ++i){
		carray[i].vext = extPotential(&carray[i]);
		carray[i].vint = pairPotential(&carray[i],&collision);
		uext += carray[i].vext;
		uint += carray[i].vint;
	}
	utot = uext + uint/2.0;
	return utot;
}

double monteCarloStep(Colloid *carray){ //returns acceptance rate
	double p=0;
	double newx = 0, newz = 0, newa = 0;
	double du = 0;	
	int collision = 0;
	int i=0;
	for(i = 0; i < N; ++i){
		newa = genrand_real2()*2.0*M_PI; //do it this way to get evenly distributed lengths
		newx = genrand_real1()*dmax; //I don't want to use another variable for this.
		newz = realZ(carray[i].z+newx*cos(newa));
		newx = realX(carray[i].x+newx*sin(newa));

		newa = carray[i].a + genrand_real2()*amax; //now the actual angle
		du = deltaU(&carray[i],newx, newz, newa, &collision);

		if( !collision && accept(du) ){
			carray[i].x = newx;
			carray[i].z = newz;
			carray[i].a = newa;
			reSortX(&carray[i]);
			reSortZ(&carray[i]);
			carray[i].vext = extPotential(&carray[i]);
			carray[i].vint = pairPotential(&carray[i],&collision);
			if(collision) printf("e");

			p += 1.0;
		}
	}
	return p/N;
}

double monteCarloSteps(Colloid *carray, int howmany){ //return acceptance rate
	double p=0;
	int i = 0;
	for(i = 0; i < howmany; ++i){
		p+=monteCarloStep(carray);
	}
	return p/howmany;
}

int accept(double du){
	if( du < 0 ) return 1;
	else{
		double paccept = exp(-du/T);
		return paccept >= genrand_real1();
	}
}

double deltaU(Colloid *c, double xnew, double znew, double anew, int* collision){
	double du = 0;
	double xtmp = 0,ztmp = 0, atmp = 0;
	xtmp = (*c).x;
	ztmp = (*c).z;
	atmp = (*c).a;
	(*c).x = xnew;
	(*c).z = znew;
	(*c).a = anew;
	reSortX(c);
	//reSortZ(c); <-- no need for this right now
	du = extPotential(c) - (*c).vext;
	du += pairPotential(c, collision) - (*c).vint;

	(*c).x = xtmp;
	(*c).z = ztmp;
	(*c).a = atmp;
	reSortX(c);
	//reSortZ(c);
	return du;
}
