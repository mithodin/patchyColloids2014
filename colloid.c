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

void insertSortedX(Colloid *list, Colloid *newitem){
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
	while( (*list).left ){
		list = (*list).left;
	}
	while( list ){
		printf("%f ",(*list).x);
		list = (*list).right;
	}
	printf("\n");
}
