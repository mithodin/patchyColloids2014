/** @file */
#include <stdlib.h>
#include <libconfig.h>
#include <string.h>
#include "config.h"
#include "load_config.h"
#include "parameters.h"

config_t *parameters; /**< parameters struct (libconfig) */
int loaded = 0; /**< count number of already loaded configurations */
config_setting_t *temperature; /**< store array of temperatures */
config_setting_t *mass2; /**< store array of masses for species 2 (twopatch) */
config_setting_t *grav; /**< store array of gravitational constants */
config_setting_t *comp; /**< store array of compositions \f$\frac{N_1}{N}\f$ */
config_setting_t *num; /**< store array of total number of particles */
int t_length = 0; /**< length of temperature */
int m2_length = 0; /**< length of mass2 */
int g_length = 0; /**< length of grav */
int comp_length = 0; /**< lenght of comp */
int num_length = 0; /**< length of num */
double x; /**< store last composition */
double lastT/** store static configurations */,lastM2/** store static configurations */,lastG/** store static configurations */,lastX/** store static configurations */,lastWidth/** store static configurations */,lastHeight/** store static configurations */,boxed; /**< store static configurations */
int lastN/** store static configurations */,lastSteps; /**< store static configurations */
const char **infile; /**< filename of initialization file */
double c_amax = 0; /**< maximum angle of rotation  (use only if loading from another simulation) */
double c_dmax = 0; /**< maximum displacement (use only if loading from another simulation) */

extern bool initdmax; /**< initialize amax and dmax? */


/**
 * Initialize the parameters struct from parameters.cfg
 * 
 * @return the parameters struct or nothing (exit the program on failure)
 */
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

/**
 * load one set of parameters (automagically traverses parameter space)
 *
 * @param c Pointer to Config Pointer. Will allocate memory there
 * @return The index of the loaded configuration (starting at 1). 0 Indicates no more parameter sets available (because of error or all sets finished)
 */
int loadParams(Config **c){ //We are relying on the fact that params is a valid config_t pointer
	int N,N1,N2,steps,iM2=0,iT=0,iG=0,iX=0,iN=0; //for intermediate storage
	double T,x,M2,g,height,width; //for intermediate storage
	const char **init=malloc(sizeof(char*)); //for reading string
	*c=(Config *)malloc(sizeof(Config)); //allocate memory for the config struct
	if(!loaded){ //The parameters have never been read
		if(config_lookup_string(parameters,"init",init) == CONFIG_FALSE){ //try to read how to initialize
			infile=NULL; //using random initialization
		}else{
			if(strcmp(*init,"boxed") == 0){
				if(config_lookup_float(parameters,"sep",&boxed) == CONFIG_FALSE){
					printf("No valid value found for sep.\n");
					return 0;
				}
				infile=NULL; //assign NULL to indicate no initialization file
			}else if(strcmp(*init,"file") == 0){
				infile=malloc(sizeof(char*)); //allocate memory for the initialization file
				if(config_lookup_string(parameters,"initfile",infile) == CONFIG_FALSE){
					printf("No valid filename given for initfile.\n");
					return 0;
				}
				boxed=0;
			}else{
				boxed=0;
				infile=NULL; //assign NULL to indicate no initalization file
			}
		}
		if(config_lookup_int(parameters,"N",&N) == CONFIG_FALSE){
			num = config_lookup(parameters,"N");
			if(num == NULL ){
				num_length = 1;
				N = -1; //-1 indicates no value found
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
				x = -1; //-1 indicates no value found
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
				printf("Invalid x. Need valid values for at least two of N,N1,N2,x. Aborting.\n");
				return 0;
			}
			if(N1+N2 == -2){
				N1 = N*x;
				N2 = N-N1;
			}else if(N+N1 == -2){
				if(x >= 1){
					printf("Invalid x for given N2. Need valid values for at least two of N,N1,N2,x. Aborting.\n");
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
				printf("Invalid N1 or N2. Need valid values for at least two of N,N1,N2,x. Aborting.\n");
				return 0;
			}else{
				N = N1+N2;
			}
		}else if(N1 == -1){
			if(N2 < 0 || N < N2){
				printf("Invalid N2 or N. Need valid values for at least two of N,N1,N2,x. Aborting.\n");
				return 0;
			}else{
				N1 = N-N2;
			}
		}else if(N2 == -1){
			if(N < N1){
				printf("Invalid N (<N1). Need valid values for at least two of N,N1,N2,x. Aborting.\n");
				return 0;
			}else{
				N2 = N-N1;
			}
		}else if(N != N1+N2){
			printf("Incorrect sum (N1+N2 != N) (N1 = %d, N2 = %d). Need valid values for at least two of N,N1,N2,x. Aborting.\n",N1,N2);
			return 0;
		}
		x = (x == -1 ? 1.0*N1/(1.0*N) : x);
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

		if(config_lookup_float(parameters,"dmax",&c_dmax) == CONFIG_FALSE){
			initdmax = true;
		}else{
			if(config_lookup_float(parameters,"amax",&c_dmax) == CONFIG_FALSE){
				initdmax = true;
				printf("Notice: no amax given, using standard initDmax\n");
			}else{
				initdmax = false;
			}
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
	(*c)->boxed = boxed;
	(*c)->dmax = c_dmax;
	(*c)->amax = c_amax;
	if( infile == NULL ){
		(*c)->loadInit = false;
	}else{
		strcpy((*c)->initIn,*infile);
		(*c)->loadInit = true;
	}
	lastT = T;
	lastM2 = M2;
	lastG = g;
	lastX = x;
	lastN = N;
	lastHeight = height;
	lastWidth = width;
	lastSteps = steps;
	printf("Config %d loaded!\n",loaded);
	printf("N = %d, N1 = %d, N2 = %d, x = %f, height = %f, width=%f, T = %f, steps = %1.0e\n",N,N1,N2,x,height,width,T,(double)steps);
	return loaded;
}
