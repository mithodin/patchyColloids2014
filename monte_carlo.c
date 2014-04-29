#include <math.h>
#include <stdio.h>
#include <stdbool.h>
#include <sys/time.h>
#include <string.h>
#include "config.h"
#include "colloid.h"
#include "statistics.h"
#include "monte_carlo.h"
#include "parameters.h"
#include "distance.h"
#include "random.h"

#define maxspan 10

const double paccept = 0.2;
const double angularPaccept = 0.44;
const double maxEnergyDeviation = 5e-3;
const double maxAccDeviation = 1e-2;
const double defaultDmax = 1e-1;
const double defaultAmax = 1e-1*2.0/3.0*M_PI;

double avg(double *, int);
void printMovie(char *, Colloid *, Config *);
void printEnergy(char *, Config *, int);
void updateUint(Colloid *, Partners *);

double extPotential(Colloid *colloid, int *collision, Config *c){
	*collision = 0;
	if( colloid->z < 0.5 || colloid->z > c->height-0.5 || colloid->x < 0.5 || colloid->x > c->width-0.5 ){
		*collision = 1;
		return 0; //invalid
	}
	return colloid->z*(colloid->sp == THREEPATCH ? M1 : c->M2)*c->g;
}

double pairPotential(Colloid *particle, int *collision, Partners *newp){
	clearPartners(newp);
	Colloid *partner = particle;
	double x = particle->x;
	double z = particle->z;
	double u = 0;
	*collision = 0;
	while( partner->below && z - partner->below->z <= sigma+delta ){
		partner = partner->below;
		if( fabs(x - partner->x) <= sigma+delta ){
			u += pairInteraction(particle,partner,collision,newp);
			if(*collision){
				return 0; //We are in an invalid state!
			}
		}
	}
	partner = particle;
	while( partner->above && partner->above->z - z <= sigma+delta ){
		partner = partner->above;
		if( fabs(x - partner->x) <= sigma+delta ){
			u += pairInteraction(particle,partner,collision,newp);
			if(*collision){
				return 0; //We are in an invalid state!
			}
		}
	}
	return u;
}

void initDmax(Colloid *carray, Config *c){
	c->dmax = defaultDmax;
	c->amax = defaultAmax;
	FILE *out = fopen(c->out,"a+");
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
		fprintf(out,"Uint: %f Avg: %f paccept: %f\n",c->Uint,avg(u,i),pnow);
		fflush(out);
	}while( ( i < maxspan || fabs(c->Uint/avg(u,i) - 1.0) > maxEnergyDeviation ) && i < 500);
	fprintf(out,"Equilibrium reached\n");
	fflush(out);

	d = c->dmax;
	c->dmax = 0.0;
	do{
		pnow = monteCarloSteps(carray,4000,c,NULL,NULL);
		c->amax *= pnow/angularPaccept;
		fprintf(out,"paccept: %f, amax = %e\n",pnow,c->amax);
		fflush(out);
	}while(fabs(pnow/angularPaccept - 1.0) > maxAccDeviation && c->amax <= 2.0*M_PI);
	fprintf(out,"Found amax = %e\n",c->amax);
	fflush(out);

	c->dmax = d;
	do{
		pnow = monteCarloSteps(carray,4000,c,NULL,NULL);
		c->dmax *= pnow/paccept;
		fprintf(out,"paccept: %f, dmax = %e\n",pnow,c->dmax);
		fflush(out);
	}while(fabs(pnow/paccept - 1.0) > maxAccDeviation);
	fprintf(out,"Using dmax = %e, amax = %e PI at paccept = %f\n",c->dmax,c->amax/M_PI,pnow);
	fflush(out);
	fclose(out);
}

double totalEnergy(Colloid *carray, Config *c){ //Give an Array here!
	double utot = 0;
	c->Uext = 0;
	c->Uint = 0;
	int collision = 0;
	int i = 0;
	for(i = 0;i < c->N; ++i){
		carray[i].vext = extPotential(&carray[i],&collision,c);
		carray[i].vint = pairPotential(&carray[i],&collision,carray[i].partners);
		c->Uext += carray[i].vext;
		c->Uint += carray[i].vint;
	}
	checkAllBonds(carray,c);
	c->Uint /= 2.0;
	utot = c->Uext + c->Uint;
	return utot;
}

void updateUint(Colloid *c, Partners *newp){
	int i;
	for(i = 0; i < patches(c); ++i){
		if( c->partners->partners[i] != NULL ){
			breakBond(c,i);
		}
	}
	for(i = 0; i < patches(c); ++i){
		if( newp->partners[i] != NULL ){
			newBond(c,newp->partners[i],i,newp->site[i]);
		}
	}
}

double monteCarloStep(Colloid *carray, Config *c, Stats *stats){ //returns acceptance rate
	double p=0;
	double oldx = 0, oldz = 0, olda = 0;
	double du = 0, duint = 0, duext = 0;	
	Partners newp;
	int collision = 0;
	int i;
	for(i = 0; i < c->N; ++i){
		oldx = carray[i].x;
		oldz = carray[i].z;
		olda = carray[i].a;

		carray[i].z = carray[i].z+c->dmax*(dsfmt_genrand_open_close(&(c->myrand))*2.0-1.0);
		carray[i].x = carray[i].x+c->dmax*(dsfmt_genrand_open_close(&(c->myrand))*2.0-1.0);
		carray[i].a = fmod(carray[i].a + (2.0*dsfmt_genrand_open_open(&(c->myrand))-1.0)*(c->amax), 2.0*M_PI);

		reSortZ(&carray[i]);

		du = deltaU(&carray[i], &duint, &duext, &newp, &collision, c);

		if( collision || !accept(du,c) ){
			carray[i].x = oldx;
			carray[i].z = oldz;
			carray[i].a = olda;
			reSortZ(&carray[i]);
		}else{
			carray[i].vext += duext;
			c->Utot += du;
			c->Uint += duint;
			c->Uext += duext;
			updateUint(&carray[i],&newp);
#ifdef PM_DEBUG
			if ( (carray[i].z - oldz != 0 && carray[i].x - oldx != 0) || carray[i].a - olda != 0 ) carray[i].haveMoved = true;
#endif
			p += 1.0;
		}
		if ( stats ){
			updateDensity(carray[i].z,carray[i].sp,c,stats);
			updateF(carray[i].z,carray[i].vint/(-U0),carray[i].sp,c,stats);
		}
	}
	if ( stats ){
		++(stats->samplingCount);
	}
	return p/(c->N);
}

double monteCarloSteps(Colloid *carray, int howmany, Config *c, Stats *stats, FILE *out){ //return acceptance rate
	bool verbose = false;
	double p=0;
	int i = 0;
	char movieFile[60];
	char energyFile[60];
	if ( stats ){
		sprintf(movieFile,"movie-%s",c->posOut);
		sprintf(energyFile,"energy-%s",c->statOut);
	}
	if(c->simRate != 0){
		double sETA = ((double)howmany)/c->simRate;
		if( out && sETA > 10){ verbose = true; }
		int hours = (int)floor(sETA/60.0/60.0);
		sETA -= hours*60*60;
		int minutes = (int)floor(sETA/60.0);
		sETA -= minutes*60;
		if ( verbose ) {
			fprintf(out,"Running %1.0e steps (eta: %dh %dmin %ds)\n[",(double)howmany,hours,minutes,(int)ceil(sETA));
			fflush(out);
		}
	}else if (out){
		verbose = true;
		fprintf(out,"Running %1.0e steps\n[",(double)howmany);
		fflush(out);
	}
	struct timeval start,stop;
	gettimeofday(&start,NULL);
	int onePerc = howmany/100;
	int k = 0, j = 0;
	for(i = 0; i < 100; ++i){
		if( verbose ){
			if( i%10 == 0){
				fprintf(out,"%d%%",i);
			}else{
				fprintf(out,".");
			}
			fflush(out);
		}
		for(k = 0; k < onePerc; ++k){
			if ( j%100 == 0 ){
				p+=monteCarloStep(carray,c,stats);
				if ( stats && j%2000 == 0 ){
					printMovie(movieFile,carray,c);
					printEnergy(energyFile,c,j);
				}
			}else{
				p+=monteCarloStep(carray,c,NULL);
			}
			++j;
		}
		checkAllBonds(carray,c);
	}
	gettimeofday(&stop,NULL);
	if (verbose ){
		fprintf(out,"100%%]\n");
		fflush(out);
	}
	double secs = (stop.tv_sec - start.tv_sec) + (double)(stop.tv_usec - start.tv_usec)/10e6;
	c->simRate=((double)howmany)/secs;
	int hours = (int)floor(secs/60.0/60.0);
	secs -= hours*60*60;
	int minutes = (int)floor(secs/60.0);
	secs -= minutes*60;
	if (verbose ){
		fprintf(out,"Time elapsed: %dh %dmin %fs\n",hours,minutes,secs);
		fflush(out);
	}
#ifdef PM_DEBUG
	if ( !out ){
		out = stdout;
	}
	int notmoved = 0;
	for(i=0;i<c->N;++i){
		if( !carray[i].haveMoved ){
			++notmoved;
			printf("vint: %lf, vext: %lf\n",carray[i].vint,carray[i].vext);
		}
		carray[i].haveMoved = false;
	}
	printf("%d particles have never been moved in %d rounds.\n",notmoved,howmany);
	fflush(stdout);
#endif
	return p/howmany;
}

int accept(double du, Config *c){
	if( du < 0 ) return 1;
	else{
		double paccept = exp(-du/c->T);
		paccept -= dsfmt_genrand_close_open(&(c->myrand));
		return paccept >= 0;
	}
}

double deltaU(Colloid *colloid, double *duint, double *duext, Partners *newp, int *collision, Config *c){
	*collision = 0;
	*duext = extPotential(colloid,collision,c) - colloid->vext;
	if(*collision) return 0; //We will reject this move anyway
	*duint = pairPotential(colloid,collision,newp) - colloid->vint;
	return *duint+*duext;
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

void printMovie(char *movieFile, Colloid *particles, Config *c){
	FILE *file = fopen(movieFile,"a");
	fprintf(file,"%d\n",c->N);
	fprintf(file,"frame\n");
	int i;
	for(i=0;i<c->N;++i){
		fprintf(file,"%s\t%f\t%f\t0\n",particles[i].sp == THREEPATCH?"C":"N",particles[i].x,particles[i].z);
	}
	fclose(file);
}

void printEnergy(char *energyFile, Config *c, int step){
	FILE *file = fopen(energyFile,"a");
	fprintf(file,"%d\t%e\t%e\t%e\n",step,c->Utot,c->Uext,c->Uint);
	fclose(file);
}
