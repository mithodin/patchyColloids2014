#include <stdlib.h>
#include <libconfig.h>
#include <string.h>
#include "load_config.h"
#include "parameters.h"

config_t *parameters;
int loaded = 0;
config_setting_t *temperature;
config_setting_t *mass2;
config_setting_t *grav;
int t_length = 0;
int m2_length = 0;
int g_length = 0;

extern char fn[40];
extern char statFn[40];
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
			mass2 = config_lookup(parameters,"M2");
			if(mass2 == NULL){
				printf("No valid value found for M2. Aborting.\n");
				return 0;
			}else{
				M2 = config_setting_get_float_elem(mass2, loaded);
				m2_length = config_setting_length(mass2);
			}
		}else{ m2_length = 1; }
		if(config_lookup_float(parameters,"g",&g) == CONFIG_FALSE){
			grav = config_lookup(parameters,"g");
			if(grav == NULL){
				printf("No valid value found for g. Aborting.\n");
				return 0;
			}else{
				g = config_setting_get_float_elem(grav, loaded);
				g_length = config_setting_length(grav);
			}
		}else{ g_length = 1; }

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
		}else{ t_length = 1; }
		sprintf(fn,"positions-T%d-M%d-G%d.dat",0,0,0);
		sprintf(statFn,"statistics-T%d-M%d-G%d.dat",0,0,0);
		loaded = 1;
		return loaded;
	}else if(loaded < t_length * m2_length * g_length){
		int iM2 = loaded%m2_length;
		int iT = (loaded/m2_length)%t_length;
		int iG = (loaded/m2_length/t_length)%g_length;
		T = temperature ? config_setting_get_float_elem(temperature, iT) : T;
		M2 = mass2 ? config_setting_get_float_elem(mass2, iM2) : M2;
		g = grav ? config_setting_get_float_elem(grav, iG) : g;
		sprintf(fn,"positions-T%d-M%d-G%d.dat",iT,iM2,iG);
		sprintf(statFn,"statistics-T%d-M%d-G%d.dat",iT,iM2,iG);
		++loaded;
		return loaded;
	}
	return 0;
}
