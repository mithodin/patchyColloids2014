typedef enum {TWOPATCH,THREEPATCH} species;

struct colloid {
	double x; //Coordinates
	double z;
	double a;
	double vext; //current external energy
	double vint; //current internal energy
	species sp; //what species?
	struct colloid *left; //Keeps a sorted double linked list in x-direction
	struct colloid *right;
	struct colloid *above; //and one in z-direction
	struct colloid *below;
};

typedef struct colloid Colloid;

void newColloid(species, Colloid *);
void insertSortedX(Colloid *list, Colloid *newitem);
void insertSortedZ(Colloid *list, Colloid *newitem);
void insertLeft(Colloid *list, Colloid *newitem);
void insertRight(Colloid *list, Colloid *newitem);
void insertBelow(Colloid *list, Colloid *newitem);
void insertAbove(Colloid *list, Colloid *newitem);
void printColloidsSortedX(Colloid *);
void printColloidsSortedZ(Colloid *);
void makePeriodicX(Colloid *);
void makePeriodicZ(Colloid *);
void swapRight(Colloid *);
void swapLeft(Colloid *);
void swapUp(Colloid *);
void swapDown(Colloid *);
void reSortX(Colloid *);
void reSortZ(Colloid *);

double colloidDistance(Colloid *, Colloid *);
double patchPositionX(Colloid *, int);
double patchPositionZ(Colloid *, int);
int patches(Colloid *);
double pairInteraction(Colloid *, Colloid *, int *);
int collisions(Colloid *);
