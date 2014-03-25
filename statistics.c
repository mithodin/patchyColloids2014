#include <stdio.h>
#include <stdlib.h>

#include "config.h"
#include "colloid.h"
#include "parameters.h"
#include "statistics.h"

int binIndex(double,double,int);
void norm(Stats *stat);

Stats *initStats(int bin){
	Stats *stat = malloc(sizeof(Stats));
	stat->bins = bin;
	stat->rho1 = (double *)calloc(bin,sizeof(double));
	stat->rho2 = (double *)calloc(bin,sizeof(double));
	stat->f1 = (double *)calloc(bin,sizeof(double));
	stat->f2 = (double *)calloc(bin,sizeof(double));
	stat->M1 = 0;
	stat->M2 = 0;
	return stat;
}

void printStats(Colloid *carray, double height, Stats *stat, char *statFn){
	norm(stat);
	FILE *statFile = fopen(statFn,"w");
	if( statFile ){
		int i;
		fprintf(statFile,"#density profile\n#Position(middle of bin)\trho1\trho2\n");
		for(i = 0; i < stat->bins; ++i){
			fprintf(statFile,"%f\t%f\t%f\n",(i+0.5)/stat->bins*height,stat->rho1[i],stat->rho2[i]);
		}
		fprintf(statFile,"\n\n");

		fprintf(statFile,"#bonds profile\n#Position(middle of bin)\tf1\tf2\n");
		for(i = 0; i < stat->bins; ++i){
			fprintf(statFile,"%f\t%f\t%f\n",(i+0.5)/stat->bins*height,stat->f1[i],stat->f2[i]);
		}
	}else{
		printf("error writing statFile - %s\n",statFn);
	}
	fclose(statFile);
}

void norm(Stats *stat){
	int i;
	for(i = 0;i<stat->bins;++i){
		if ( stat->rho1[i] != 0 ) stat->f1[i] /= 3.0*stat->rho1[i];
		if ( stat->rho2[i] != 0 ) stat->f2[i] /= 2.0*stat->rho2[i];
		stat->rho1[i] /= stat->M1;
		stat->rho2[i] /= stat->M2;
	}
}

int binIndex(double z, double height, int bins){
	int i = (int)(z/height*bins);
	i = i < bins ? i : bins-1;
	return i;
}

void updateDensity(double z, species kind, Config *c, Stats *stat){
	int i = binIndex(z,c->height,stat->bins);
	switch(kind){
		case THREEPATCH:
			stat->rho1[i] += 1.0;
			++(stat->M1);
			break;
		case TWOPATCH:
			stat->rho2[i] += 1.0;
			++(stat->M2);
			break;
	}
}

void updateF(double z, double val, species kind, Config *c, Stats *stat){
	int i = binIndex(z,c->height,stat->bins);
	switch(kind){
		case THREEPATCH:
			stat->f1[i] += val;
			break;
		case TWOPATCH:
			stat->f2[i] += val;
			break;
	}
}
