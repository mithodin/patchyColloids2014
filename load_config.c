#include <stdlib.h>
#include <libconfig.h>
#include <string.h>
#include "config.h"
#include "load_config.h"
#include "parameters.h"

config_t *parameters;
int loaded = 0;
config_setting_t *temperature;
config_setting_t *mass2;
config_setting_t *grav;
config_setting_t *comp; //Composition
config_setting_t *num;
int t_length = 0;
int m2_length = 0;
int g_length = 0;
int comp_length = 0;
int num_length = 0;
double x;
double lastT,lastM2,lastG,lastX,lastWidth,lastHeight;
int lastN,lastSteps;


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

int loadParams(Config **c){ //We are relying on the fact that params is a valid config_t pointer
	int N,N1,N2,steps,iM2=0,iT=0,iG=0,iX=0,iN=0;
	double T,x,M2,g,height,width;
	*c=(Config *)malloc(sizeof(Config));
	if(!loaded){
		if(config_lookup_int(parameters,"N",&N) == CONFIG_FALSE){
			num = config_lookup(parameters,"N");
			if(num == NULL ){
				num_length = 1;
				N = -1;
			}else{
				N = config_setting_get_int_elem(num, loaded);
				num_length = config_setting_length(num);
			}
		}else{ num_length = 1; }
		if(config_lookup_int(parameters,"N1",&N1) == CONFIG_FALSE) N1 = -1;
		if(config_lookup_int(parameters,"N2",&N2) == CONFIG_FALSE) N2 = -1;
		if(config_lookup_float(parameters,"x",&x) == CONFIG_FALSE){
			comp = config_lookup(parameters,"x");
			if(comp == NULL){
				comp_length = 1;
				x = -1;
			}else{
				x = config_setting_get_float_elem(comp, loaded);
				comp_length = config_setting_length(comp);
			}
		}else{ comp_length = 1; }

		if( N1+N2+N == -3 || N1+N2+x == -3 || N1+N+x == -3 || N2+N+x == -3 ){
			printf("Need valid values for at least two of N,N1,N2,x. Aborting.\n");
			return 0;
		}
		if(x != -1){
			if(x < 0 || x > 1){
				printf("1 Need valid values for at least two of N,N1,N2,x. Aborting.\n");
				return 0;
			}
			if(N1+N2 == -2){
				N1 = N*x;
				N2 = N*(1-x);
			}else if(N+N1 == -2){
				if(x >= 1){
					printf("2 Need valid values for at least two of N,N1,N2,x. Aborting.\n");
					return 0;
				}else{
					N1 = (1.0/(1.0-x)-1.0)*N2;
				}
			}else if(N+N2 == -2){
				N2 = (1.0/x-1.0)*N1;
			}
		}
		if(N == -1){
			if(N1 < 0 || N2 < 0){
				printf("3 Need valid values for at least two of N,N1,N2,x. Aborting.\n");
				return 0;
			}else{
				N = N1+N2;
			}
		}else if(N1 == -1){
			if(N2 < 0 || N < N2){
				printf("4 Need valid values for at least two of N,N1,N2,x. Aborting.\n");
				return 0;
			}else{
				N1 = N-N2;
			}
		}else if(N2 == -1){
			if(N < N1){
				printf("5 Need valid values for at least two of N,N1,N2,x. Aborting.\n");
				return 0;
			}else{
				N2 = N-N1;
			}
		}else if(N != N1+N2){
			printf("x = %f, N=%d, N1=%d, N2=%d\n",x,N,N1,N2);
			printf("6 Need valid values for at least two of N,N1,N2,x. Aborting.\n");
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
		loaded = 1;
	}else if(loaded < t_length * m2_length * g_length * comp_length * num_length){
		iM2 = loaded%m2_length;
		iT = (loaded/m2_length)%t_length;
		iG = (loaded/m2_length/t_length)%g_length;
		iX = (loaded/m2_length/t_length/g_length)%comp_length;
		iN = (loaded/m2_length/t_length/g_length/comp_length)%num_length;
		T = temperature ? config_setting_get_float_elem(temperature, iT) : lastT;
		M2 = mass2 ? config_setting_get_float_elem(mass2, iM2) : lastM2;
		g = grav ? config_setting_get_float_elem(grav, iG) : lastG;
		x = comp ? config_setting_get_float_elem(comp, iX) : lastX;
		N = num ? config_setting_get_int_elem(num, iN) : lastN;
		if( x < 0 || x > 1 ){
			printf("Invalid x value!\n");
			return 0;
		}
		N1 = N*x;
		N2 = N*(1-x);
		width = lastWidth;
		height = lastHeight;
		steps = lastSteps;
		++loaded;
	}else{
		return 0;
	}
	sprintf((*c)->posOut,"positions-T%d-M%d-G%d-X%d-N%d.dat",iT,iM2,iG,iX,iN);
	sprintf((*c)->statOut,"statistics-T%d-M%d-G%d-X%d-N%d.dat",iT,iM2,iG,iX,iN);
	(*c)->N = N;
	(*c)->N1 = N1;
	(*c)->N2 = N2;
	(*c)->M2 = M2;
	(*c)->height = height;
	(*c)->width = width;
	(*c)->T = T;
	(*c)->steps = steps;
	(*c)->g = g;
	(*c)->simRate = 0;
	lastT = T;
	lastM2 = M2;
	lastG = g;
	lastX = x;
	lastN = N;
	lastHeight = height;
	lastWidth = width;
	lastSteps = steps;
	printf("Config %d loaded!\n",loaded);
	return loaded;
}
