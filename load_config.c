#include <stdlib.h>
#include <libconfig.h>
#include "headers.h"

config_t *getParams(void){ //loads the params from a file
	config_t *parameters=malloc(sizeof(struct config_t));

	FILE *paramFile = fopen("parameters.cfg","r");
	if(paramFile == NULL){
		return NULL;
	}

	config_init(parameters);
	if(config_read(parameters, paramFile)) return parameters;
	return NULL;
}

//Maybe I will extend this to handle arrays of parameters.
void loadParams(config_t *params){ //We are relying on the fact that params is a valid config_t pointer
	if(config_lookup_int(params,"N",&N) == CONFIG_FALSE) N = -1;
	if(config_lookup_int(params,"N1",&N1) == CONFIG_FALSE) N1 = -1;
	if(config_lookup_int(params,"N2",&N2) == CONFIG_FALSE) N2 = -1;

	if(N == -1){
		if(N1 == -1 || N2 == -1){
			printf("Need valid values for at least two of N,N1,N2. Aborting.\n");
			exit(1);
		}else{
			N = N1+N2;
		}
	}else if(N1 == -1){
		if(N2 == -1 || N < N2){
			printf("Need valid values for at least two of N,N1,N2. Aborting.\n");
			exit(1);
		}else{
			N1 = N-N2;
		}
	}else if(N2 == -1){
		if(N < N1){
			printf("Need valid values for at least two of N,N1,N2. Aborting.\n");
			exit(1);
		}else{
			N2 = N-N1;
		}
	}else if(N != N1+N2){
		printf("Need valid values for at least two of N,N1,N2. Aborting.\n");
		exit(1);
	}
	if(config_lookup_float(params,"U0",&U0) == CONFIG_FALSE){
		printf("No valid value found for U0. Aborting.\n");
		exit(1); //We need this value!
	}
	if(config_lookup_float(params,"M2",&M2) == CONFIG_FALSE){
		printf("No valid value found for M2. Aborting.\n");
		exit(1);
	}
}
