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
	col->partners->partners[0] = NULL;
	col->partners->partners[1] = NULL;
	col->partners->partners[2] = NULL;
	col->haveMoved = false;
}

void checkBonds(Colloid *c){
	int i;
	for(i=0;i<patches(c);++i){
		if(c->partners->partners[i]){
			if(c->partners->partners[i]->partners->partners[c->partners->site[i]] != c){
				printf("unreciprocated bond!\n");
				exit(-1);
			}
		}
	}
}

void checkAllBonds(Colloid *particles, Config *c){
	int i;
	for(i=0;i<c->N;++i){
		checkBonds(&particles[i]);
	}
}

void clearPartners(Partners *p){
	int i;
	for(i=0;i<3;++i){
		p->partners[i]=NULL;
		p->site[i]=-1;
	}
}

void newBond(Colloid *c1, Colloid *c2, int site1, int site2){
	if(c1->partners->partners[site1] || c2->partners->partners[site2]){
		printf("error. bond sites not clear!\n");
		exit(-1);
	}
	c1->partners->partners[site1] = c2;
	c1->partners->site[site1] = site2;
	c2->partners->partners[site2] = c1;
	c2->partners->site[site2] = site1;
	c2->vint -= U0;
	c1->vint -= U0;
}

void breakBond(Colloid *c1, int site1){
	c1->partners->partners[site1]->partners->partners[c1->partners->site[site1]] = NULL;
	c1->partners->partners[site1]->partners->site[c1->partners->site[site1]] = -1;
	c1->partners->partners[site1]->vint += U0;

	c1->partners->partners[site1] = NULL;
	c1->partners->site[site1] = -1;
	c1->vint += U0;
}

//For z direction
void insertSortedZ(Colloid *list, Colloid *newitem){
	while( list->below && list->below->z > newitem->z ){
		list=list->below;
	}
	while( list->above && list->above->z < newitem->z ){
		list=list->above;
	}
	if( list->z > newitem->z ){
		insertBelow(list, newitem);
	}else{
		insertAbove(list, newitem);
	}
}

void insertAbove(Colloid *list, Colloid *newitem){
	Colloid *tmp = list->above;
	list->above = newitem;
	newitem->below = list;
	newitem->above = tmp;
	if( tmp ){ //If tmp is not a NULL pointer
		tmp->below = newitem;
	}
}

void insertBelow(Colloid *list, Colloid *newitem){
	Colloid *tmp = list->below;
	list->below = newitem;
	newitem->above = list;
	newitem->below = tmp;
	if( tmp ){
		tmp->above = newitem;
	}
}

void printColloidsSortedZ(Colloid *list){
	printf("Colloids sorted by z coordinate:\n");
	while( list->below ){
		list = list->below;
	}
	int count = 0;
	do {
		if( list->sp == THREEPATCH ){
			printf("(%f) ",list->z);
		}else{
			printf("[%f] ",list->z);
		}
		list = list->above;
		++count;
	}while( list );
	printf("\nTotal Number: %d\n",count);
}

void reSortZ(Colloid *list, Config *c){ //Point to the changed element!
	while( list->above && list->z > list->above->z ){
		swapUp(list);
	}
	while( list->below && list->z < list->below->z ){
		swapDown(list);
	}
}

void swapUp(Colloid *list){ //assumes there is a particle above
	Colloid *curAbove = list->above;
	Colloid *curBelow = list->below;
	list->above = curAbove->above;
	list->below = curAbove;
	if( curBelow ) curBelow->above = curAbove;
	if( curAbove->above ) curAbove->above->below = list;
	curAbove->below = curBelow;
	curAbove->above = list;
}

void swapDown(Colloid *list){ //assumes there is a particle below
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
			for( j = 0; j < patches(c2); ++j ){
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
	return distance(c1->x,c1->z,c2->x,c2->z);
}

double patchPositionX(Colloid *c, int which){
	if( (*c).sp == THREEPATCH ){
		return (*c).x + sigma/2.0*cos((*c).a + (which%3)*2.0/3.0*M_PI);
	}else if( (*c).sp == TWOPATCH ){
		return (*c).x + sigma/2.0*cos((*c).a + (which%2)*M_PI);
	}else{
		printf("invalid patch position.\n");
		exit(-1);
	}
}

double patchPositionZ(Colloid *c, int which){
	if( (*c).sp == THREEPATCH ){
		return (*c).z + sigma/2.0*sin((*c).a + (which%3)*2.0/3.0*M_PI);
	}else if( (*c).sp == TWOPATCH ){
		return (*c).z + sigma/2.0*sin((*c).a + (which%2)*M_PI);
	}else{
		printf("invalid patch position.\n");
		exit(-1);
	}
}

int patches(Colloid *c){
	switch( c->sp ){
		case THREEPATCH: return 3; break;
		case TWOPATCH: return 2; break;
		default:
			printf("error\n");
			exit(-1);
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
