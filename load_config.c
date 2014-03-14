#include <stdlib.h>
#include <libconfig.h>
#include <string.h>
#include "load_config.h"
#include "parameters.h"

config_t *parameters;
int loaded = 0;
config_setting_t *temperature;
int t_length = 0;

extern char fn[40];
extern FILE *initFile;

config_t *getParams(void){ //loads the params from a file
	parameters=(struct config_t*)malloc(sizeof(struct config_t));

	FILE *paramFile = fopen("parameters.cfg","r");
	if(paramFile == NULL){
		printf("Error: No file 'parameters.cfg' in 'bin' directory\n");
		exit(1);
	}

	config_init(parameters);
	if(config_read(parameters, paramFile)) return parameters;
	else{
		printf("Error loading parameters file. Check syntax.\n");
		exit(1);
	}
}

//Maybe I will extend this to handle arrays of parameters.
int loadParams(){ //We are relying on the fact that params is a valid config_t pointer
	if(!loaded){
		const char *infile;
		if(config_lookup_string(parameters,"init",&infile) == CONFIG_TRUE){
			if(strcmp(infile,"random") != 0){
				initFile = fopen(infile, "r");
			}
		}
		if(config_lookup_int(parameters,"N",&N) == CONFIG_FALSE) N = -1;
		if(config_lookup_int(parameters,"N1",&N1) == CONFIG_FALSE) N1 = -1;
		if(config_lookup_int(parameters,"N2",&N2) == CONFIG_FALSE) N2 = -1;

		if(N == -1){
			if(N1 == -1 || N2 == -1){
				printf("Need valid values for at least two of N,N1,N2. Aborting.\n");
				return 0;
			}else{
				N = N1+N2;
			}
		}else if(N1 == -1){
			if(N2 == -1 || N < N2){
				printf("Need valid values for at least two of N,N1,N2. Aborting.\n");
				return 0;
			}else{
				N1 = N-N2;
			}
		}else if(N2 == -1){
			if(N < N1){
				printf("Need valid values for at least two of N,N1,N2. Aborting.\n");
				return 0;
			}else{
				N2 = N-N1;
			}
		}else if(N != N1+N2){
			printf("Need valid values for at least two of N,N1,N2. Aborting.\n");
			return 0;
		}
		if(config_lookup_float(parameters,"M2",&M2) == CONFIG_FALSE){
			printf("No valid value found for M2. Aborting.\n");
			return 0;
		}
		if(config_lookup_float(parameters,"height",&height) == CONFIG_FALSE){
			printf("No valid value found for height. Aborting.\n");
			return 0;
		}
		if(config_lookup_float(parameters,"width",&width) == CONFIG_FALSE){
			printf("No valid value found for width. Aborting.\n");
			return 0;
		}
		if(config_lookup_int(parameters,"steps",&steps) == CONFIG_FALSE) steps = 1000;
		if(config_lookup_float(parameters,"T",&T) == CONFIG_FALSE){
			temperature = config_lookup(parameters,"T");
			if(temperature == NULL){
				printf("No valid value found for T. Aborting.\n");
				return 0;
			}else{
				T = config_setting_get_float_elem(temperature, loaded);
				t_length = config_setting_length(temperature);
			}
		}
		sprintf(fn,"positions-T%d.dat",loaded);
		loaded = 1;
		return 1;
	}else if(loaded < t_length){
		T = config_setting_get_float_elem(temperature, loaded);
		sprintf(fn,"positions-T%d.dat",loaded);
		++loaded;
		return 1;
	}
	return 0;
}
