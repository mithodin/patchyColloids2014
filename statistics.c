#include <stdio.h>
#include <stdlib.h>

#include "colloid.h"
#include "parameters.h"
#include "statistics.h"

char statFn[40];
double *rho1;
double *rho2;
double *f1;
double *f2;
int bins;

int M1dens = 0, M2dens = 0, M1f = 0, M2f = 0;

int binIndex(double);
void norm(void);

void initStats(int bin){
	bins = bin;
	rho1 = (double *)calloc(bins,sizeof(double));
	rho2 = (double *)calloc(bins,sizeof(double));
	f1 = (double *)calloc(bins,sizeof(double));
	f2 = (double *)calloc(bins,sizeof(double));
}

void printStats(Colloid *carray){
	norm();
	FILE *stats = fopen(statFn,"w");
	if( stats ){
		int i;
		fprintf(stats,"#density profile\n#Position(middle of bin)\trho1\trho2\n");
		for(i = 0; i < bins; ++i){
			fprintf(stats,"%f\t%f\t%f\n",(i+0.5)/bins*height,rho1[i],rho2[i]);
		}
		fprintf(stats,"\n\n");

		fprintf(stats,"#bonds profile\n#Position(middle of bin)\tf1\tf2\n");
		for(i = 0; i < bins; ++i){
			fprintf(stats,"%f\t%f\t%f\n",(i+0.5)/bins*height,f1[i],f2[i]);
		}
	}else{
		printf("error writing stats - %s\n",statFn);
	}
	fclose(stats);
}

void norm(void){
	int i;
	for(i = 0;i<bins;++i){
		rho1[i] /= M1dens;
		rho2[i] /= M2dens;
		if ( rho1[i] != 0 ) f1[i] /= 3.0*rho1[i]*M1f;
		if ( rho2[i] != 0 ) f2[i] /= 2.0*rho2[i]*M2f;
	}
}

int binIndex(double z){
	int i = (int)(z/height*bins);
	i = i < bins ? i : bins-1;
	return i;
}

void updateDensity(double z, species kind){
	int i = binIndex(z);
	switch(kind){
		case THREEPATCH:
			rho1[i] += 1.0;
			++M1dens;
			break;
		case TWOPATCH:
			rho2[i] += 1.0;
			++M2dens;
			break;
	}
}

void updateF(double z, double val, species kind){
	int i = binIndex(z);
	switch(kind){
		case THREEPATCH:
			f1[i] += val;
			++M1f;
			break;
		case TWOPATCH:
			f2[i] += val;
			++M2f;
			break;
	}
}
