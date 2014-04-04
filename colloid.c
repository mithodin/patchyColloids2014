#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "config.h"
#include "colloid.h"
#include "parameters.h"
#include "distance.h"

void newColloid(species sp, Colloid *col){
	col->above = NULL;
	col->below = NULL;
	col->vext = 0;
	col->vint = 0;
	col->sp = sp;
	col->partners = malloc(sizeof(Partners));
}

void newBond(Colloid *c1, Colloid *c2, int site1, int site2){
	c1->partners->partners[site1] = c2;
	c1->partners->site[site1] = site2;
	c2->partners->partners[site2] = c1;
	c1->partners->site[site2] = site1;
}

void breakBond(Colloid *c1, Colloid *c2, int site1, int site2){
	c1->partners->partners[site1] = NULL;
	c2->partners->partners[site2] = NULL;
}

//For z direction
void insertSortedZ(Colloid *list, Colloid *newitem){ //DO NOT USE WHEN DLL HAS BEEN MADE PERIODIC!
	while( (*list).below && (*(*list).below).z > (*newitem).z ){
		list=(*list).below;
	}
	while( (*list).above && (*(*list).above).z < (*newitem).z ){
		list=(*list).above;
	}
	if( (*list).z > (*newitem).z ){
		insertBelow(list, newitem);
	}else{
		insertAbove(list, newitem);
	}
}

void insertAbove(Colloid *list, Colloid *newitem){
	Colloid *tmp = (*list).above;
	(*list).above = newitem;
	(*newitem).below = list;
	(*newitem).above = tmp;
	if( tmp ){ //If tmp is not a NULL pointer
		(*tmp).below = newitem;
	}
}

void insertBelow(Colloid *list, Colloid *newitem){
	Colloid *tmp = (*list).below;
	(*list).below = newitem;
	(*newitem).above = list;
	(*newitem).below = tmp;
	if( tmp ){
		(*tmp).above = newitem;
	}
}

void printColloidsSortedZ(Colloid *list){
	printf("Colloids sorted by z coordinate:\n");
	while( (*list).below ){
		list = (*list).below;
	}
	do {
		if( (*list).sp == THREEPATCH ){
			printf("(%f) ",(*list).z);
		}else{
			printf("[%f] ",(*list).z);
		}
		list = (*list).above;
	}while( list );
	printf("\n");
}

void reSortZ(Colloid *list, Config *c){ //Point to the changed element!
	while( list->above && list->z > list->above->z ){
		swapUp(list);
	}
	while( list->below && list->z < list->below->z ){
		swapDown(list);
	}
}

void swapUp(Colloid *list){
	Colloid *curAbove = (*list).above;
	Colloid *curBelow = (*list).below;
	(*list).above = (*curAbove).above;
	(*list).below = curAbove;
	if( curBelow ) (*curBelow).above = curAbove;
	if( curAbove->above ) (*(*curAbove).above).below = list;
	(*curAbove).below = curBelow;
	(*curAbove).above = list;
}

void swapDown(Colloid *list){
	Colloid *curAbove = (*list).above;
	Colloid *curBelow = (*list).below;
	(*list).below = (*curBelow).below;
	(*list).above = curBelow;
	if( curBelow->below ) (*(*curBelow).below).above = list;
	if( curAbove ) (*curAbove).below = curBelow;
	(*curBelow).above = curAbove;
	(*curBelow).below = list;
}

double pairInteraction(Colloid *c1, Colloid *c2, int *collision, Partners *newp){
	double d = colloidDistance(c1,c2);
	*collision = 0;
	if(d > sigma+delta){ return 0.0; }
	else if(d < sigma ){ 
		*collision = 1;
		return 0.0;
	}
	else{
		int i = 0;
		int j = 0;
		for(i = 0; i < patches(c1); i++){
			for( j = 0; j < patches(c2); j++ ){
				double p1X = patchPositionX(c1,i);
				double p1Z = patchPositionZ(c1,i);
				double p2X = patchPositionX(c2,j);
				double p2Z = patchPositionZ(c2,j);
				
				d = distance(p1X,p1Z,p2X,p2Z);
				if( d <= delta ){
					newp->partners[i]=c2;
					newp->site[i]=j;
					return -U0;
				}
			}
		}
		return 0;
	}
}

double colloidDistance(Colloid *c1, Colloid *c2){
	return distance((*c1).x,(*c1).z,(*c2).x,(*c2).z);
}

double patchPositionX(Colloid *c, int which){
	if( (*c).sp == THREEPATCH ){
		return (*c).x + sigma/2.0*cos((*c).a + (which%3)*2.0/3.0*M_PI);
	}else if( (*c).sp == TWOPATCH ){
		return (*c).x + sigma/2.0*cos((*c).a + (which%2)*M_PI);
	}else{
		return -1; //ERROR!
	}
}

double patchPositionZ(Colloid *c, int which){
	if( (*c).sp == THREEPATCH ){
		return (*c).z + sigma/2.0*sin((*c).a + (which%3)*2.0/3.0*M_PI);
	}else if( (*c).sp == TWOPATCH ){
		return (*c).z + sigma/2.0*sin((*c).a + (which%2)*M_PI);
	}else{
		return -1; //ERROR!
	}
}

int patches(Colloid *c){
	switch( (*c).sp ){
		case THREEPATCH: return 3; break;
		case TWOPATCH: return 2; break;
		default: return 0;
	}
}

int collisions(Colloid *carray, Config *c){
	int i=0,j=0;
	for(i = 0; i<c->N; ++i){
		for(j = i+1; j < c->N; j++){
			if(colloidDistance(&carray[i], &carray[j]) < sigma ){
				return 1;
			}
		}
	}
	return 0;
}
