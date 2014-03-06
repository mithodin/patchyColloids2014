void init_genrand(unsigned long s);
void init_by_array(unsigned long init_key[], int key_length);
unsigned long genrand_int32(void);
long genrand_int31(void);
double genrand_real1(void);
double genrand_real2(void);
double genrand_real3(void);
double genrand_real4(void);
double genrand_res53(void);
unsigned long int random_seed(void);
void box_muller(double *n1, double *n2);

void init(void);
int overlap(int i);
void writepos(void);
