double extPotential(Colloid *);
double pairPotential(Colloid *, int*);
double totalEnergy(Colloid *);
double monteCarloStep(Colloid *);
double monteCarloSteps(Colloid *, int);
double deltaU(Colloid*, double, double, double, int*);
int accept(double);
