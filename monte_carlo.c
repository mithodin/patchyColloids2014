#include <math.h>
#include <stdio.h>
#include "colloid.h"
#include "monte_carlo.h"
#include "parameters.h"
#include "distance.h"
#include "mt19937ar.h"

const double paccept = 0.6;
double dmax = 5; //completely random numbers
double amax = 2.0/3.0*M_PI;

double extPotential(Colloid *c){
	return 0;
}

double pairPotential(Colloid *particle, int *collision){
	Colloid *partner = particle;
	double x = (*particle).x;
	double z = (*particle).z;
	double u = 0;
	int col = 0;
	*collision = 0;
	while( fabs(realDx(x - (*(*partner).left).x)) <= sigma+delta ){
		partner = (*partner).left;
		if( fabs(realDz(z - (*partner).z)) <= sigma+delta ){
			u += pairInteraction(particle,partner,&col);
			if(col){
			 	*collision = 1;
				return 0; //We are in an invalid state!
			}
		}
	}
	partner = particle;
	while( fabs(realDx(x - (*(*partner).right).x)) <= sigma+delta ){
		partner = (*partner).right;
		if( fabs(realDz(z - (*partner).z)) <= sigma+delta ){
			u += pairInteraction(particle,partner,&col);
			if(col){
				*collision = 1;
				return 0; //We are in an invalid state!
			}
		}
	}
	return u;
}

void initDmax(Colloid *carray){
	double pnow = monteCarloSteps(carray,1000);
	/*
	do{
		dmax = pnow/paccept*dmax;
		printf("trying dmax = %f...",dmax);
		fflush(stdout);
		pnow = monteCarloSteps(carray,1000);
		printf("paccept: %f\n",pnow);
	}while(fabs(pnow/paccept - 1.0) > 10e-2); 
*/
	printf("Using dmax = %f, amax = %f PI at paccept = %f\n",dmax,amax/M_PI,pnow);
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
	double r = 0, a = 0, oldx = 0, oldz = 0, olda = 0;
	double du = 0;	
	int collision = 0;
	int i=0;
	for(i = 0; i < N; ++i){
		oldx = carray[i].x;
		oldz = carray[i].z;
		olda = carray[i].a;

		a = genrand_real2()*2.0*M_PI; //do it this way to get evenly distributed lengths
		r = genrand_real1()*dmax;

		carray[i].z = realZ(carray[i].z+r*cos(a));
		carray[i].x = realX(carray[i].x+r*sin(a));
		carray[i].a = carray[i].a + genrand_real2()*amax;

		reSortX(&carray[i]);
		reSortZ(&carray[i]);

		du = deltaU(&carray[i], &collision);

		if( collision || !accept(du) ){
			carray[i].x = oldx;
			carray[i].z = oldz;
			carray[i].a = olda;
			reSortX(&carray[i]);
			reSortZ(&carray[i]);
		}else{
			carray[i].vext = extPotential(&carray[i]);
			carray[i].vint = pairPotential(&carray[i],&collision);
			p += 1.0;
		}
	}
	return p/N;
}

double monteCarloSteps(Colloid *carray, int howmany){ //return acceptance rate
	FILE *output = fopen("/dev/null","w");
	if(howmany > 1000){
		output = stdout;
	}
	fprintf(output,"Running %d steps: [",howmany);
	double p=0;
	int i = 0;
	for(i = 0; i < howmany; ++i){
		p+=monteCarloStep(carray);
		if(fmod((100.0*i)/howmany,1) == 0){
			if(fmod((100.0*i)/howmany,10) == 0){
				fprintf(output,"%d%%",(int)((100.0*i)/howmany));
			}else{
				fprintf(output,".");
			}
			fflush(output);
		}
	}
	fprintf(output,"100%%]\n");
	return p/howmany;
}

int accept(double du){
	if( du < 0 ) return 1;
	else{
		double paccept = exp(-du/T);
		return paccept >= genrand_real1();
	}
}

double deltaU(Colloid *c, int *collision){
	double du = 0;
	du = extPotential(c) - (*c).vext;
	du += pairPotential(c, collision) - (*c).vint;

	return du;
}
