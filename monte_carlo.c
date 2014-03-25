#include <math.h>
#include <stdio.h>
#include <stdbool.h>
#include <sys/time.h>
#include "config.h"
#include "colloid.h"
#include "statistics.h"
#include "monte_carlo.h"
#include "parameters.h"
#include "distance.h"
#include "mt19937ar.h"

const double paccept = 0.2;
const double angularPaccept = 0.4;
const double maxEnergyDeviation = 5e-3;
const double maxAccDeviation = 1e-2;
const int maxspan = 10;

double avg(double *, int);

double extPotential(Colloid *colloid, int *collision, Config *c){
	*collision = 0;
	if( colloid->z < 0.5 || colloid->z > c->height-0.5 ){
		*collision = 1;
		return 0; //invalid
	}
	return colloid->z*(colloid->sp == THREEPATCH ? M1 : c->M2)*c->g;
}

double pairPotential(Colloid *particle, int *collision, Config *c){
	Colloid *partner = particle;
	double x = (*particle).x;
	double z = (*particle).z;
	double u = 0;
	int col = 0;
	*collision = 0;
	while( fabs(realDx(x - (*(*partner).left).x, c->width)) <= sigma+delta ){
		partner = (*partner).left;
		if( fabs(z - (*partner).z) <= sigma+delta ){
			u += pairInteraction(particle,partner,&col,c);
			if(col){
			 	*collision = 1;
				return 0; //We are in an invalid state!
			}
		}
	}
	partner = particle;
	while( fabs(realDx(x - (*(*partner).right).x, c->width)) <= sigma+delta ){
		partner = (*partner).right;
		if( fabs(z - (*partner).z) <= sigma+delta ){
			u += pairInteraction(particle,partner,&col,c);
			if(col){
				*collision = 1;
				return 0; //We are in an invalid state!
			}
		}
	}
	return u;
}

void initDmax(Colloid *carray, Config *c, FILE *out){
	double d = c->dmax; //just storing this for later reference
	double u[maxspan] = {0};
	double pnow = paccept;
	int i=-1;
	do{
		++i;
		d = c->dmax;
		c->dmax = 0;
		pnow = monteCarloSteps(carray,4000,c,NULL,NULL);
		c->amax *= pnow/paccept;
		c->amax = c->amax > 2.0*M_PI ? 2.0*M_PI : c->amax; //for high temperatures, amax is practically meaningless
		c->dmax = d;
		pnow = monteCarloSteps(carray,4000,c,NULL,NULL);
		c->dmax *= pnow/paccept;
		c->Utot = totalEnergy(carray,c);
		u[i%maxspan] = c->Uint;
		//fprintf(out,"Uint: %f Avg: %f paccept: %f\n",c->Uint,avg(u,i),pnow);
	}while( ( i < maxspan || fabs(c->Uint/avg(u,i) - 1.0) > maxEnergyDeviation ) && i < 500);
	//fprintf(out,"Equilibrium reached\n");

	d = c->dmax;
	c->dmax = 0.0;
	do{
		pnow = monteCarloSteps(carray,4000,c,NULL,NULL);
		c->amax *= pnow/angularPaccept;
		//fprintf(out,"paccept: %f, amax = %e\n",pnow,c->amax);
	}while(fabs(pnow/angularPaccept - 1.0) > maxAccDeviation && c->amax <= 2.0*M_PI);
	//fprintf(out,"Found amax = %e\n",c->amax);

	c->dmax = d;
	do{
		pnow = monteCarloSteps(carray,4000,c,NULL,NULL);
		c->dmax *= pnow/paccept;
		//fprintf(out,"paccept: %f, dmax = %e\n",pnow,c->dmax);
	}while(fabs(pnow/paccept - 1.0) > maxAccDeviation);
	//fprintf(out,"Using dmax = %e, amax = %e PI at paccept = %f\n",c->dmax,c->amax/M_PI,pnow);
}

double totalEnergy(Colloid *carray, Config *c){ //Give an Array here!
	double utot = 0;
	c->Uext = 0;
	c->Uint = 0;
	int collision = 0;
	int i = 0;
	for(i = 0;i < c->N; ++i){
		carray[i].vext = extPotential(&carray[i],&collision,c);
		carray[i].vint = pairPotential(&carray[i],&collision,c);
		c->Uext += carray[i].vext;
		c->Uint += carray[i].vint;
	}
	c->Uint /= 2.0;
	utot = c->Uext + c->Uint;
	return utot;
}

double monteCarloStep(Colloid *carray, Config *c, Stats *stats){ //returns acceptance rate
	double p=0;
	double oldx = 0, oldz = 0, olda = 0;
	double du = 0;	
	int collision = 0;
	int i=0;
	for(i = 0; i < c->N; ++i){
		oldx = carray[i].x;
		oldz = carray[i].z;
		olda = carray[i].a;

		carray[i].z = realZ(carray[i].z+c->dmax*(genrand_real1()*2.0-1.0),c->height);
		carray[i].x = realZ(carray[i].x+c->dmax*(genrand_real1()*2.0-1.0),c->width);
		carray[i].a = carray[i].a + (2.0*genrand_real2()-1.0)*(c->amax);

		reSortX(&carray[i], c);
		reSortZ(&carray[i], c);

		du = deltaU(&carray[i], &collision, c);

		if( collision || !accept(du,c->T) ){
			carray[i].x = oldx;
			carray[i].z = oldz;
			carray[i].a = olda;
			reSortX(&carray[i],c);
			reSortZ(&carray[i],c);
		}else{
			carray[i].vext = extPotential(&carray[i],&collision,c);
			carray[i].vint = pairPotential(&carray[i],&collision,c);
			p += 1.0;
		}
		if ( stats ){
			updateDensity(carray[i].z,carray[i].sp,c,stats);
			updateF(carray[i].z,carray[i].vint/(-U0),carray[i].sp,c,stats);
		}
	}
	return p/c->N;
}

double monteCarloSteps(Colloid *carray, int howmany, Config *c, Stats *stats, FILE *out){ //return acceptance rate
	bool verbose = false;
	double p=0;
	int i = 0;
	if(c->simRate != 0){
		double sETA = ((double)howmany)/c->simRate;
		if(sETA > 10){ verbose = true; }
		int hours = (int)floor(sETA/60.0/60.0);
		sETA -= hours*60*60;
		int minutes = (int)floor(sETA/60.0);
		sETA -= minutes*60;
		//if ( verbose ) //fprintf(out,"Running %1.0e steps (eta: %dh %dmin %ds)\n[",(double)howmany,hours,minutes,(int)ceil(sETA));
	}else{
		verbose = true;
		//fprintf(out,"Running %1.0e steps\n[",(double)howmany);
	}
	struct timeval start,stop;
	gettimeofday(&start,NULL);
	int onePerc = howmany/100;
	int k = 0, j = 0;
	for(i = 0; i < 100; ++i){
		if( verbose ){
			if( i%10 == 0){
				//fprintf(out,"%d%%",i);
			}else{
				//fprintf(out,".");
			}
			//fflush(out);
		}
		for(k = 0; k < onePerc; ++k){
			if ( j%100 == 0 ){
				p+=monteCarloStep(carray,c,stats);
			}else{
				p+=monteCarloStep(carray,c,NULL);
			}
			++j;
		}
	}
	gettimeofday(&stop,NULL);
	//if (verbose ) //fprintf(out,"100%%]\n");
	double secs = (stop.tv_sec - start.tv_sec) + (double)(stop.tv_usec - start.tv_usec)/10e6;
	c->simRate=((double)howmany)/secs;
	int hours = (int)floor(secs/60.0/60.0);
	secs -= hours*60*60;
	int minutes = (int)floor(secs/60.0);
	secs -= minutes*60;
	//if (verbose ) //fprintf(out,"Time elapsed: %dh %dmin %fs\n",hours,minutes,secs);
	return p/howmany;
}

int accept(double du, double T){
	if( du < 0 ) return 1;
	else{
		double paccept = exp(-du/T);
		return paccept >= genrand_real1();
	}
}

double deltaU(Colloid *colloid, int *collision, Config *c){
	double du = 0;
	int col = 0;
	*collision = 0;
	du = extPotential(colloid,&col,c) - colloid->vext;
	*collision += col;
	du += pairPotential(colloid, &col,c) - colloid->vint;
	*collision += col;
	return du;
}

double avg(double *array, int index){
	int span = (index+1)>maxspan ? maxspan : index+1;
	int i = 0;
	double avg = 0;
	for(;i < span; ++i){
		avg += array[i];
	}
	return avg/span;
}
