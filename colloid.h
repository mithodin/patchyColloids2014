typedef enum {TWOPATCH,THREEPATCH} species;

typedef struct {
	double x; //Coordinates
	double z;
	double a;
	species sp; //what species?
} Colloid;
