struct stat{
	long double *rho1;
	long double *rho2;
	long double *f1;
	long double *f2;
	int bins;
	long samplingCount;
};

typedef struct stat Stats;

void printStats(Colloid *, Stats*, Config *);
Stats *initStats(int);
void updateDensity(double, species, Config*, Stats*);
void updateF(double, double, species, Config*, Stats*);
