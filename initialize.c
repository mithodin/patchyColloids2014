#include <stdlib.h>
#include <math.h>
#include "mt19937ar.h"
#include "parameters.h"
#include "colloid.h"
#include "distance.h"
#include "initialize.h"

void initParticles(Colloid *particles){
	//particles=malloc(sizeof(Colloid)*N);
	int i=0;
	Colloid *tmp = NULL;
	Colloid *this = NULL;
	while(i<N){
		particles[i].x = genrand_real2()*width;
		particles[i].z = genrand_real2()*height;
		particles[i].a = genrand_real2()*2*M_PI;
		if(noCollision(i,particles[i].x,particles[i].z,particles)){
			this = &particles[i];
			if(i<N1){
				newColloid(THREEPATCH,this);
			}else{
				newColloid(TWOPATCH,this);
			}
			if(i>0){
				insertSortedX(&particles[i-1], this);
			}
			 ++i;
		}
	}
}

int noCollision(int i, double x, double z, Colloid *particles){
	int j=0;
	for(j=0;j<N;j++){
		if(j==i) continue;
		if(distance(x,z,particles[j].x,particles[j].z) <= 1.0) return 0;
	}
	return 1;
}
