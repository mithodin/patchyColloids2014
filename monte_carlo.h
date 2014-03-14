double extPotential(Colloid *);
double pairPotential(Colloid *, int*);
double totalEnergy(Colloid *, double *, double *);
double monteCarloStep(Colloid *);
double monteCarloSteps(Colloid *, int);
double deltaU(Colloid*, int*);
int accept(double);
void initDmax(Colloid *);
