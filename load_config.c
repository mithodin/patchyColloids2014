#include "headers.h"

config_t *getParams(void){
	config_t *parameters=malloc(sizeof(struct config_t));

	FILE *paramFile = fopen("parameters.cfg","r");
	if(paramFile == NULL){
		return NULL;
	}

	config_init(parameters);
	if(config_read(parameters, paramFile)) return parameters;
	return NULL;
}
