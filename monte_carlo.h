/**
 * @file monte_carlo.h
 * @brief Defines all the functions for the Monte Carlo simulation
 */
double extPotential(Colloid *, int*, Config *);
double pairPotential(Colloid *, int*, Partners *);
double totalEnergy(Colloid *, Config*);
double monteCarloStep(Colloid *, Config*, Stats *);
double monteCarloSteps(Colloid *, int, Config*, Stats *, FILE*);
double deltaU(Colloid*, double *, double *, Partners *, int*, Config*);
int accept(double, Config*);
void initDmax(Colloid*, Config*);
