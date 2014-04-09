#import <stdbool.h>

typedef enum {TWOPATCH,THREEPATCH} species;

struct partners {
	struct colloid *partners[3];
	int site[3];
};

typedef struct partners Partners;

struct colloid {
	double x; //Coordinates
	double z;
	double a;
	double vext; //current external energy
	double vint; //current internal energy
	species sp; //what species?
	struct colloid *above; //and one in z-direction
	struct colloid *below;
	Partners *partners; //all bonding partners. 2 is always to be NULL for TWOPATCH
	bool haveMoved;
};

typedef struct colloid Colloid;

void checkBonds(Colloid *);
void checkAllBonds(Colloid *, Config *);
void newColloid(species, Colloid *);
void clearPartners(Partners *);
void newBond(Colloid *, Colloid *, int, int);
void breakBond(Colloid *, int);
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
double pairInteraction(Colloid *, Colloid *, int *, Partners *);
int collisions(Colloid *, Config *c);
