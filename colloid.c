#include <stdlib.h>
#include <stdio.h>
#include "colloid.h"

void newColloid(species sp, Colloid *col){
	(*col).left = NULL;
	(*col).right = NULL;
	(*col).above = NULL;
	(*col).below = NULL;
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
	while( (*list).below && (*(*list).below).z < (*list).z ){
		list = (*list).below;
	}
	do {
		if( (*list).sp == THREEPATCH ){
			printf("(%f) ",(*list).x);
		}else{
			printf("[%f] ",(*list).x);
		}
		list = (*list).above;
	}while( list && (*(*list).below).z < (*list).z );
	printf("\n");
}

void makePeriodicZ(Colloid *list){
	Colloid *first=list;
	Colloid *last=list;
	while( (*first).below ){
		first = (*first).below;
	}
	while( (*last).above ){
		last = (*last).above;
	}
	(*first).below = last;
	(*last).above = first;
}

void reSortZ(Colloid *list){ //Point to the changed element!
	while( (*list).z > (*(*list).above).z && (*(*list).above).z > (*(*list).below).z ){
		swapUp(list);
	}
	while( (*list).z < (*(*list).below).z && (*(*list).above).z > (*(*list).below).z ){
		swapDown(list);
	}
}

void swapUp(Colloid *list){
	Colloid *curAbove = (*list).above;
	Colloid *curBelow = (*list).below;
	(*list).above = (*curAbove).above;
	(*list).below = curAbove;
	(*curBelow).above = curAbove;
	(*(*curAbove).above).below = list;
	(*curAbove).below = curBelow;
	(*curAbove).above = list;
}
