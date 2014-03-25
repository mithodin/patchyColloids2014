double extPotential(Colloid *, int*, double);
double pairPotential(Colloid *, int*);
double totalEnergy(Colloid *, double *, double *, Config*);
double monteCarloStep(Colloid *, Config*, Stats *);
double monteCarloSteps(Colloid *, int, Config*, Stats *, FILE*);
double deltaU(Colloid*, int*, Config*);
int accept(double, double);
void initDmax(Colloid*, Config*, FILE*);
