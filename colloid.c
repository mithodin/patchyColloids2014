/** @file */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "config.h"
#include "colloid.h"
#include "parameters.h"
#include "distance.h"

/**
 * Initialize a new Colloid struct
 *
 * @param sp Define the species of the new colloid
 * @param col Pointer to Colloid struct. Memory must be assigned beforehand!
 */
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

/**
 * Check if all bonds of the colloid are reciprocal
 *
 * Will exit the program with nonzero status if invalid bonds exist
 *
 * @param c Colloid to check
 */
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

/**
 * Check all bonds of a colloid array
 *
 * @param particles Colloid array
 * @param c Configuration struct
 */
void checkAllBonds(Colloid *particles, Config *c){
	int i;
	for(i=0;i<c->N;++i){
		checkBonds(&particles[i]);
	}
}

/**
 * Set all bonds to NULL
 *
 * @param p Pointer to Partners struct. Pass colloid->partners here
 */
void clearPartners(Partners *p){
	int i;
	for(i=0;i<3;++i){
		p->partners[i]=NULL;
		p->site[i]=-1;
	}
}

/**
 * Creates a new bond between two colloids
 *
 * Will exit the program if one of the bonding sites is not clear
 * Bonding site definitions are not being checked!
 *
 * @param c1 Colloid 1
 * @param c2 Colloid 2
 * @param site1 Bonding site for first colloid
 * @param site2 Bonding site for second colloid
 */
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

/**
 * Break an existing bond
 *
 * @param c1 One of the bonding partners
 * @param site1 Bonding site on given colloid. Checks are not done that this is valid
 */
void breakBond(Colloid *c1, int site1){
	c1->partners->partners[site1]->partners->partners[c1->partners->site[site1]] = NULL;
	c1->partners->partners[site1]->partners->site[c1->partners->site[site1]] = -1;
	c1->partners->partners[site1]->vint += U0;

	c1->partners->partners[site1] = NULL;
	c1->partners->site[site1] = -1;
	c1->vint += U0;
}

/**
 * Insert a colloid into a dl list sorted by z coordinate
 *
 * @param list Entry point into the doubly linked list. May be any colloid in the list
 * @param newitem Colloid to add to the list
 */
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

/**
 * Insert new colloid above given entry point
 *
 * @param newitem Colloid to add
 * @param list newitem will be included above this item
 */
void insertAbove(Colloid *list, Colloid *newitem){
	Colloid *tmp = list->above;
	list->above = newitem;
	newitem->below = list;
	newitem->above = tmp;
	if( tmp ){ //If tmp is not a NULL pointer
		tmp->below = newitem;
	}
}

/**
 * Insert new colloid below given entry point
 *
 * @param newitem Colloid to add
 * @param list newitem will be included below this item
 */
void insertBelow(Colloid *list, Colloid *newitem){
	Colloid *tmp = list->below;
	list->below = newitem;
	newitem->above = list;
	newitem->below = tmp;
	if( tmp ){
		tmp->above = newitem;
	}
}

/**
 * Print a list of colloids sorted by z coordinate
 *
 * @param list Entry point to the list. May be anywhere in the list
 */
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

/**
 * Re-Sort a dl list of colloids
 *
 * @param list Must point to the changed colloid!
 */ 
void reSortZ(Colloid *list){ //Point to the changed element!
	while( list->above && list->z > list->above->z ){
		swapUp(list);
	}
	while( list->below && list->z < list->below->z ){
		swapDown(list);
	}
}

/**
 * Swap the colloid one up
 */
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

/**
 * Swap the colloid one down
 */
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

/**
 * Calculate bonding energy between two colloids
 *
 * Will not update c1 or c2 in any way. Detected bonds go to newp.
 *
 * @param c1 Colloid 1
 * @param c2 Colloid 2
 * @param collision Pointer to collision indicator (1 = collision, 0 = no collision)
 * @param newp Pointer to bond storage
 * @return Returns the bonding energy
 */
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

/**
 * Calculate distance between two colloids (center to center)
 *
 * @param c1 Colloid 1
 * @param c2 Colloid 2
 * @return Returns the distance between the centers
 */
double colloidDistance(Colloid *c1, Colloid *c2){
	return distance(c1->x,c1->z,c2->x,c2->z);
}

/**
 * Calculate the x coordinate of a patch
 *
 * If patch index is invalid, this will exit the program
 *
 * @param c Colloid struct
 * @param which Which patch? (0,1,2 for threepatch, 0,1 for twopatch)
 * @return The x coordinate of the patch
 */
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

/**
 * Calculate the z coordinate of a patch
 *
 * If patch index is invalid, this will exit the program
 *
 * @param c Colloid struct
 * @param which Which patch? (0,1,2 for threepatch, 0,1 for twopatch)
 * @return The z coordinate of the patch
 */
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

/**
 * Number of patches on a colloid
 *
 * @param c Colloid
 * @return Number of patches on the colloid
 */
int patches(Colloid *c){
	switch( c->sp ){
		case THREEPATCH: return 3; break;
		case TWOPATCH: return 2; break;
		default:
			printf("error\n");
			exit(-1);
	}
}

/**
 * Detect any collisions between colloids
 *
 * @param carray Array of colloids of length c->N
 * @param c Configuration
 * @return 0 for no collision, 1 for 1 ore more collisions
 */
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
