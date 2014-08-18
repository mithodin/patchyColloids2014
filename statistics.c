/** @file */
#include <stdio.h>
#include <stdlib.h>

#include "config.h"
#include "colloid.h"
#include "parameters.h"
#include "statistics.h"

int binIndex(double,double,int);
void norm(Stats *stat, Config *c);

/**
 * Create a new Stats struct
 *
 * @param bin Number of bins
 * @return Pointer to the newly created Stats struct
 */
Stats *initStats(int bin){
	Stats *stat = malloc(sizeof(Stats));
	stat->bins = bin;
	stat->rho1 = (long double *)calloc(bin,sizeof(long double));
	stat->rho2 = (long double *)calloc(bin,sizeof(long double));
	stat->f1 = (long double *)calloc(bin,sizeof(long double));
	stat->f2 = (long double *)calloc(bin,sizeof(long double));
	stat->samplingCount = 0;
	return stat;
}

/**
 * Print the collected statistics to c->statOut
 *
 * @param stat Pointer to Stats struct. Normalize before printing!
 * @param c Configuration struct
 */
void printStats(Stats *stat, Config *c){
	norm(stat, c);
	FILE *statFile = fopen(c->statOut,"w");
	if( statFile ){
		int i;
		fprintf(statFile,"#density profile\n#Position(middle of bin)\trho1\trho2\n");
		for(i = 0; i < stat->bins; ++i){
			fprintf(statFile,"%f\t%Lf\t%Lf\n",(i+0.5)/stat->bins*c->height,stat->rho1[i],stat->rho2[i]);
		}
		fprintf(statFile,"\n\n");

		fprintf(statFile,"#bonds profile\n#Position(middle of bin)\tf1\tf2\n");
		for(i = 0; i < stat->bins; ++i){
			fprintf(statFile,"%f\t%Lf\t%Lf\n",(i+0.5)/stat->bins*c->height,stat->f1[i],stat->f2[i]);
		}
	}else{
		printf("error writing statFile - %s\n",c->statOut);
	}
	fclose(statFile);
}

/**
 * Normalize the collected statistics.
 * Do this only after the simulation is over!
 *
 * @param stat Statistics struct
 * @param c Configuration struct
 */
void norm(Stats *stat, Config *c){
	int i;
	for(i = 0;i<stat->bins;++i){
		if ( stat->rho1[i] != 0 ) stat->f1[i] /= 3.0*stat->rho1[i];
		if ( stat->rho2[i] != 0 ) stat->f2[i] /= 2.0*stat->rho2[i];
		stat->rho1[i] /= stat->samplingCount*c->width*(c->height/stat->bins);
		stat->rho2[i] /= stat->samplingCount*c->width*(c->height/stat->bins);
	}
}

/**
 * Calculate the correct bin for a z coordinate
 *
 * @param z z Coordinate
 * @param height Height of the simulation box
 * @param bins Number of bins
 * @return Index of correct bin
 */
int binIndex(double z, double height, int bins){
	int i = (int)(z/height*bins);
	i = i < bins ? i : bins-1;
	return i;
}

/**
 * Update the density array
 *
 * @param z z Coordinate of current particle
 * @param kind Which species does the particle belong to
 * @param c Configuration struct
 * @param stat Statistics struct
 */
void updateDensity(double z, species kind, Config *c, Stats *stat){
	int i = binIndex(z,c->height,stat->bins);
	switch(kind){
		case THREEPATCH:
			stat->rho1[i] += 1.0;
			break;
		case TWOPATCH:
			stat->rho2[i] += 1.0;
			break;
	}
}

/**
 * Update the bonds array
 *
 * @param z z Coordinate of current particle
 * @param val How many bonds does the particle have
 * @param kind Which species does the particle belong to
 * @param c Configuration struct
 * @param stat Statistics struct
 */
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
