/** @file
 * @brief This file defines the Colloid struct and the functions to handle colloids
 */
#import <stdbool.h>

/**
 * Enum to indicate a colloid's species
 *
 * The values should be self-explanatory
 */
typedef enum {TWOPATCH,THREEPATCH} species;

struct partners {
	struct colloid *partners[3]; /**< Array of bonded colloids. Index indicates bonding site */
	int site[3]; /**< Array of partner's bonding sites */
};

typedef struct partners Partners;

struct colloid {
	double x; /**< x coordinate */
	double z; /**< z coordinate */
	double a; /**< angle in radians */
	double vext; /**< current external energy */
	double vint; /**< current internal energy */
	species sp; /**< species of the colloid (twopatch/threepatch) */
	struct colloid *above; /**< which colloid is next above this one */
	struct colloid *below; /**< which colloid is next below this one */
	Partners *partners; /**< bonding partners. See Partners */
	bool haveMoved; /**< true if this colloid has ever made a successfull MC step */
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
void swapUp(Colloid *);
void swapDown(Colloid *);
void reSortZ(Colloid *);

double colloidDistance(Colloid *, Colloid *);
double patchPositionX(Colloid *, int);
double patchPositionZ(Colloid *, int);
int patches(Colloid *);
double pairInteraction(Colloid *, Colloid *, int *, Partners *);
int collisions(Colloid *, Config *c);
