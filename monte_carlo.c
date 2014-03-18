#include <math.h>
#include <stdio.h>
#include <sys/time.h>
#include "colloid.h"
#include "monte_carlo.h"
#include "parameters.h"
#include "distance.h"
#include "mt19937ar.h"

const double paccept = 0.2;
const double angularPaccept = 0.4;
const double maxEnergyDeviation = 5e-3;
const double maxAccDeviation = 1e-2;
double dmax = 1.0; //completely random numbers
double amax = 0.1*2.0/3.0*M_PI;
double simRate = 0; //mc steps per second.
const int maxspan = 10;

double avg(double *, int);

double extPotential(Colloid *c, int *collision){
	*collision = 0;
	if( c->z < 0.5 || c->z > height-0.5 ){
		*collision = 1;
		return 0; //invalid
	}
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
		if( fabs(z - (*partner).z) <= sigma+delta ){
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
		if( fabs(z - (*partner).z) <= sigma+delta ){
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
	double d = dmax; //just storing this for later reference
	double u[300] = {0}; //I assume this will be enough
	double pnow = paccept;
	int i=-1;
	do{
		++i;
		d = dmax;
		dmax = 0;
		pnow = monteCarloSteps(carray,4000);
		amax *= pnow/paccept;
		dmax = d;
		pnow = monteCarloSteps(carray,4000);
		dmax *= pnow/paccept;
		Utot = totalEnergy(carray,&Uext,&Uint);
		u[i] = Uint;
		printf("Uint: %f Avg: %f paccept: %f\n",Uint,avg(u,i),pnow);
	}while(i < maxspan || fabs(Uint/avg(u,i) - 1.0) > maxEnergyDeviation);
	printf("Equilibrium reached\n");

	d = dmax;
	dmax = 0.0;
	do{
		pnow = monteCarloSteps(carray,4000);
		amax *= pnow/angularPaccept;
		printf("paccept: %f, amax = %e\n",pnow,amax);
	}while(fabs(pnow/angularPaccept - 1.0) > maxAccDeviation);
	printf("Found amax = %e\n",amax);

	dmax = d;
	do{
		pnow = monteCarloSteps(carray,4000);
		dmax *= pnow/paccept;
		printf("paccept: %f, dmax = %e\n",pnow,dmax);
	}while(fabs(pnow/paccept - 1.0) > maxAccDeviation);
	printf("Using dmax = %e, amax = %e PI at paccept = %f\n",dmax,amax/M_PI,pnow);
}

double totalEnergy(Colloid *carray, double *uext, double *uint){ //Give an Array here!
	double utot = 0;
	*uext = 0;
	*uint = 0;
	int collision = 0;
	int i = 0;
	for(i = 0;i < N; ++i){
		carray[i].vext = extPotential(&carray[i],&collision);
		carray[i].vint = pairPotential(&carray[i],&collision);
		*uext += carray[i].vext;
		*uint += carray[i].vint;
	}
	*uint /= 2.0;
	utot = *uext + *uint;
	return utot;
}

double monteCarloStep(Colloid *carray){ //returns acceptance rate
	double p=0;
	double oldx = 0, oldz = 0, olda = 0;
	double du = 0;	
	int collision = 0;
	int i=0;
	for(i = 0; i < N; ++i){
		oldx = carray[i].x;
		oldz = carray[i].z;
		olda = carray[i].a;

		carray[i].z = realZ(carray[i].z+dmax*(genrand_real1()*2.0-1.0));
		carray[i].x = realZ(carray[i].x+dmax*(genrand_real1()*2.0-1.0));
		carray[i].a = carray[i].a + (2.0*genrand_real2()-1.0)*amax;

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
			carray[i].vext = extPotential(&carray[i],&collision);
			carray[i].vint = pairPotential(&carray[i],&collision);
			p += 1.0;
		}
	}
	return p/N;
}

double monteCarloSteps(Colloid *carray, int howmany){ //return acceptance rate
	FILE *output = fopen("/dev/null","w");
	double p=0;
	int i = 0;
	if(simRate != 0){
		double sETA = ((double)howmany)/simRate;
		if(sETA > 10){ output = stdout; }
		int hours = (int)floor(sETA/60.0/60.0);
		sETA -= hours*60*60;
		int minutes = (int)floor(sETA/60.0);
		sETA -= minutes*60;
		fprintf(output, "Running %1.0e steps (eta: %dh %dmin %ds)\n [",(double)howmany,hours,minutes,(int)ceil(sETA));
	}else{
		output = stdout;
		fprintf(output,"Running %1.0e steps\n[",(double)howmany);
	}
	struct timeval start,stop;
	gettimeofday(&start,NULL);
	int onePerc = howmany/100;
	int k = 0;
	for(i = 0; i < 100; ++i){
		if(i%10 == 0){
			fprintf(output,"%d%%",i);
		}else{
			fprintf(output,".");
		}
		fflush(output);
		for(k = 0; k < onePerc; ++k){
			p+=monteCarloStep(carray);
		}
	}
	gettimeofday(&stop,NULL);
	fprintf(output,"100%%]\n");
	double secs = (stop.tv_sec - start.tv_sec) + (double)(stop.tv_usec - start.tv_usec)/10e6;
	simRate=((double)howmany)/secs;
	int hours = (int)floor(secs/60.0/60.0);
	secs -= hours*60*60;
	int minutes = (int)floor(secs/60.0);
	secs -= minutes*60;
	fprintf(output,"Time elapsed: %dh %dmin %fs\n",hours,minutes,secs);
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
	int col = 0;
	*collision = 0;
	du = extPotential(c,&col) - (*c).vext;
	*collision += col;
	du += pairPotential(c, &col) - (*c).vint;
	*collision += col;
	return du;
}

double avg(double *array, int index){
	int span = (index+1)>maxspan?maxspan:(index+1);
	int i = index - span +1;
	double avg = 0;
	for(;i <= index; ++i){
		avg += array[i];
	}
	return avg/span;
}
