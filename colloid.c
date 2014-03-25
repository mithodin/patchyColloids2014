#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "config.h"
#include "colloid.h"
#include "parameters.h"
#include "distance.h"

void newColloid(species sp, Colloid *col){
	(*col).left = NULL;
	(*col).right = NULL;
	(*col).above = NULL;
	(*col).below = NULL;
	(*col).vext = 0;
	(*col).vint = 0;
	(*col).sp = sp;
}

//For x direction
void insertSortedX(Colloid *list, Colloid *newitem){ //DO NOT USE WHEN DLL HAS BEEN MADE PERIODIC!
	while( (*list).left && (*(*list).left).x > (*newitem).x ){
		list=(*list).left;
	}
	while( (*list).right && (*(*list).right).x < (*newitem).x ){
		list=(*list).right;
	}
	if( (*list).x > (*newitem).x ){
		insertLeft(list, newitem);
	}else{
		insertRight(list, newitem);
	}
}

void insertRight(Colloid *list, Colloid *newitem){
	Colloid *tmp = (*list).right;
	(*list).right = newitem;
	(*newitem).left = list;
	(*newitem).right = tmp;
	if( tmp ){ //If tmp is not a NULL pointer
		(*tmp).left = newitem;
	}
}

void insertLeft(Colloid *list, Colloid *newitem){
	Colloid *tmp = (*list).left;
	(*list).left = newitem;
	(*newitem).right = list;
	(*newitem).left = tmp;
	if( tmp ){
		(*tmp).right = newitem;
	}
}

void printColloidsSortedX(Colloid *list){
	printf("Colloids sorted by x coordinate:\n");
	while( (*list).left && (*(*list).left).x < (*list).x ){
		list = (*list).left;
	}
	do {
		if( (*list).sp == THREEPATCH ){
			printf("(%f) ",(*list).x);
		}else{
			printf("[%f] ",(*list).x);
		}
		list = (*list).right;
	}while( list && (*(*list).left).x < (*list).x );
	printf("\n");
}

void makePeriodicX(Colloid *list){
	Colloid *first=list;
	Colloid *last=list;
	while( (*first).left ){
		first = (*first).left;
	}
	while( (*last).right ){
		last = (*last).right;
	}
	(*first).left = last;
	(*last).right = first;
}

void reSortX(Colloid *list, Config *c){ //Point to the changed element!
	while( realDx( (*list).x - (*(*list).right).x, c->width ) > 0 ){
		swapRight(list);
	}
	while( realDx( (*list).x - (*(*list).left).x, c->width ) < 0 ){
		swapLeft(list);
	}
}

void swapRight(Colloid *list){
	Colloid *curRight = (*list).right;
	Colloid *curLeft = (*list).left;
	(*list).right = (*curRight).right;
	(*list).left = curRight;
	(*curLeft).right = curRight;
	(*(*curRight).right).left = list;
	(*curRight).left = curLeft;
	(*curRight).right = list;
}

void swapLeft(Colloid *list){
	Colloid *curRight = (*list).right;
	Colloid *curLeft = (*list).left;
	(*list).left = (*curLeft).left;
	(*list).right = curLeft;
	(*(*curLeft).left).right = list;
	(*curRight).left = curLeft;
	(*curLeft).right = curRight;
	(*curLeft).left = list;
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

double pairInteraction(Colloid *c1, Colloid *c2, int *collision, Config *c){
	double d = colloidDistance(c1,c2,c);
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
				
				d = distance(p1X,p1Z,p2X,p2Z,c);
				if( d <= delta ){
					return -U0;
				}
			}
		}
		return 0;
	}
}

double colloidDistance(Colloid *c1, Colloid *c2, Config *c){
	return distance((*c1).x,(*c1).z,(*c2).x,(*c2).z,c);
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
			if(colloidDistance(&carray[i], &carray[j],c) < sigma ){
				return 1;
			}
		}
	}
	return 0;
}
