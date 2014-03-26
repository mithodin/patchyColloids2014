typedef enum {TWOPATCH,THREEPATCH} species;

struct colloid {
	double x; //Coordinates
	double z;
	double a;
	double vext; //current external energy
	double vint; //current internal energy
	species sp; //what species?
	struct colloid *above; //and one in z-direction
	struct colloid *below;
};

typedef struct colloid Colloid;

void newColloid(species, Colloid *);
void insertSortedZ(Colloid *list, Colloid *newitem);
void insertBelow(Colloid *list, Colloid *newitem);
void insertAbove(Colloid *list, Colloid *newitem);
void printColloidsSortedZ(Colloid *);
void makePeriodicZ(Colloid *);
void swapUp(Colloid *);
void swapDown(Colloid *);
void reSortZ(Colloid *, Config *c);

double colloidDistance(Colloid *, Colloid *);
double patchPositionX(Colloid *, int);
double patchPositionZ(Colloid *, int);
int patches(Colloid *);
double pairInteraction(Colloid *, Colloid *, int *);
int collisions(Colloid *, Config *c);
