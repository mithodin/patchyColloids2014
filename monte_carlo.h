#include <stdbool.h>

double extPotential(Colloid *, int*);
double pairPotential(Colloid *, int*);
double totalEnergy(Colloid *, double *, double *);
double monteCarloStep(Colloid *, bool, Stats *);
double monteCarloSteps(Colloid *, int, Stats *);
double deltaU(Colloid*, int*);
int accept(double);
void initDmax(Colloid *);
