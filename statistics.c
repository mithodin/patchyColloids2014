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

double **rho(Colloid*);
double **f(Colloid*);
int binIndex(double);

void printStats(Colloid *carray){
	bins = 10;
	double **rs = rho(carray);
	rho1 = rs[0];
	rho2 = rs[1];

	rs = f(carray);
	f1 = rs[0];
	f2 = rs[1];

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

double **rho(Colloid *carray){
	double *r1 = calloc(bins,sizeof(double));
	double *r2 = calloc(bins,sizeof(double));
	static double *rs[2];
	rs[0] = r1;
	rs[1] = r2;
	int i = 0;
	for(i = 0;i<N;++i){
		switch(carray[i].sp){
			case THREEPATCH:
				++(r1[binIndex(carray[i].z)]);
				break;
			case TWOPATCH:
				++(r2[binIndex(carray[i].z)]);
				break;
		}
	}
	for(i = 0;i<bins;++i){
		r1[i] /= N1;
		r2[i] /= N2;
	}
	return rs;
}

double **f(Colloid *carray){ //rho needs to be calculated first
	double *f1 = calloc(bins,sizeof(double));
	double *f2 = calloc(bins,sizeof(double));
	static double *fs[2];
	fs[0] = f1;
	fs[1] = f2;
	int i = 0;
	for(i = 0;i<N;++i){
		switch(carray[i].sp){
			case THREEPATCH:
				f1[binIndex(carray[i].z)] += carray[i].vint/(-U0);
				break;
			case TWOPATCH:
				f2[binIndex(carray[i].z)] += carray[i].vint/(-U0);
				break;
		}
	}
	for(i = 0;i<bins;++i){
		f1[i] /= N1*3.0*rho1[i];
		f2[i] /= N2*2.0*rho2[i];
	}
	return fs;
}

int binIndex(double z){
	int i = (int)(z/height*bins);
	i = i < bins ? i : bins-1;
	return i;
}
