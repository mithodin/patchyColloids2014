struct conf{
	*int done;
	int N;
	int N1;
	int N2;
	double M2;
	double height;
	double width;
	double T;
	int steps;
	double g;
	double dmax;
	double amax;

	char out[40];
	char statOut[40];
	char posOut[40];
};

typedef struct conf Config;
