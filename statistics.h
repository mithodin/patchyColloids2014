struct stat{
	double *rho1;
	double *rho2;
	double *f1;
	double *f2;
	int bins;
	int M1,M2;
};

typedef struct stat Stats;

void printStats(Colloid *, double, Stats*, char*);
Stats *initStats(int);
void updateDensity(double, species, Config*, Stats*);
void updateF(double, double, species, Config*, Stats*);
