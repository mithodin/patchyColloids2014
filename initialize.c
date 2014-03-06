#include <libconfig.h>
#include <stdlib.h>
#include "mt19937ar.h"
#include "parameters.h"
#include "colloid.h"

void initParticles(Colloid *particles){
	particles=malloc(sizeof(Colloid)*N);
}
