/*************************************************
 * Main Simulation file for patchy colloids
 * in a gravitational field
 *
 * Bachelor's Thesis by Lucas Treffenstädt
 *************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <libconfig.h> 	//We are using this to load our parameters from a file.

#include "mt19937ar.h" 	//Generating random numbers
#include "headers.h" 	//Global headers file

config_t *parameters;

int main(void){ //This is only for testing so far.
	parameters=getParams();
	if(parameters==NULL) return 1;
	
	char *greeting=malloc(sizeof(char)*100);
	config_lookup_string(parameters,"greeting",&greeting);
	printf(greeting);
	printf("\n");
	return 0;
}
