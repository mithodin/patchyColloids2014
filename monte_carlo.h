double extPotential(Colloid *, int*, Config *);
double pairPotential(Colloid *, int*);
double totalEnergy(Colloid *, Config*);
double monteCarloStep(Colloid *, Config*, Stats *);
double monteCarloSteps(Colloid *, int, Config*, Stats *, FILE*);
double deltaU(Colloid*, int*, Config*);
int accept(double, double);
void initDmax(Colloid*, Config*);
