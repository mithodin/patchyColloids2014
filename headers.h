// Global headers file

// For load_config.c
config_t *getParams(void); //load the configuration from the standard conf file
void loadParams(config_t *);

//All the parameters
extern int N; //N1 + N2 = N
extern int N1;
extern int N2;
extern double U0;
extern const double M1; //m1 is always 1
extern double M2; //in units of m1
