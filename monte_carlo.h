double extPotential(Colloid *, int*, Config *);
double pairPotential(Colloid *, int*, Partners *);
double totalEnergy(Colloid *, Config*);
double monteCarloStep(Colloid *, Config*, Stats *);
double monteCarloSteps(Colloid *, int, Config*, Stats *, FILE*);
double deltaU(Colloid*, double *, double *, Partners *, int*, Config*);
int accept(double, double);
void initDmax(Colloid*, Config*);
