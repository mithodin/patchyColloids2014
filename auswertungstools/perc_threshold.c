#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>

//MAKE SURE HEIGHT IS A MULTIPLE OF BINHEIGHT!
#define BINHEIGHT 2
#define HEIGHT 200
#define NUM_BINS (HEIGHT/BINHEIGHT)

typedef enum{READ_NUMBER,READ_PARTICLE,SKIP_ONE,ERROR,END_OF_FRAME} file_state;

//STATUS VARIABLES
int N; //Number of particles
int frame = 0; //How many frames have been read?
int current_particle;
double *x=NULL;
double *z=NULL;
double *patch_x[3]={NULL};
double *patch_z[3]={NULL};
unsigned short int *patches=NULL;
unsigned int *bin=NULL;
long double bin_bonds[NUM_BINS]={0.0};
unsigned long int bin_patches[NUM_BINS]={0};
long double bin_composition[NUM_BINS]={0.0};
unsigned long int bin_n[NUM_BINS]={0.0};
const double sigma = 1.0;
const double delta = 0.11965683746373795115;
//----------------

file_state read_number(const char *);
file_state read_particle(const char *);
unsigned int binof(double);
unsigned int bondsof(void);
double patch_distance(int p1, int patch1, int p2, int patch2);
double distance(int p1, int p2);
bool hasbond(int p1, int p2);

int main(int argc, char **argv){
	if(argc != 2){
		printf("usage: pbond <moviefile>\n");
		return 1;
	}

	FILE *datei=fopen(argv[1],"r");
	if(datei == NULL){
		printf("file %s does not exist.\n",argv[1]);
		return 1;
	}
	int maxlength = 100;
	char zeile[maxlength];
	file_state fs = READ_NUMBER;
	while(fgets(zeile,maxlength,datei) != NULL){
		switch(fs){
			case END_OF_FRAME:
				printf("Frame %d done.\n",++frame);
			case READ_NUMBER:
				fs = read_number(zeile);
				break;
			case SKIP_ONE:
				fs = READ_PARTICLE;
				break;
			case READ_PARTICLE:
				fs = read_particle(zeile);
				break;
			default:
				printf("error: file not valid!\n");
				return 1;
		}
	}
	fclose(datei);
	if(fs == END_OF_FRAME){
		printf("All frames (%d) done.\n",++frame);
		int i;
		for(i=0;i<NUM_BINS;++i){
			if(bin_patches[i]>0){
				bin_bonds[i]/=bin_patches[i];
			}else{
				printf("%lu patches in bin %d\n",bin_patches[i],i);
			}
			if(bin_n[i]>0){
				bin_composition[i]/=bin_n[i];
			}else{
				printf("%lu particles in bin %d\n",bin_n[i],i);
			}
		}
		printf("writing..."); fflush(stdout);
		FILE *output=fopen("pbonds.dat","w");
		fprintf(output,"#bin index\tz(bin)\tpbond\tcomposition\n");
		for(i=0;i<NUM_BINS;++i){
			fprintf(output,"%d\t%lf\t%Lf\t%Lf\n",i,(double)((i+0.5)*BINHEIGHT),bin_bonds[i],bin_composition[i]);
		}
		fflush(output);
		fclose(output);
		printf("done.\n");
		free(x);
		free(z);
		free(patches);
		free(bin);
		for(i=0;i<3;++i){
			free(patch_x[i]);
			free(patch_z[i]);
		}
		return 0;
	}else{
		return 1;
	}
}

file_state read_number(const char *zeile){
	int vars=0;
	if((vars=sscanf(zeile,"%d",&N))!=1){
		printf("Read %d vars instead of 1 (N).\n%s\n",vars,zeile);
		return ERROR;
	}
	int i;
	free(x);
	free(z);
	free(patches);
	free(bin);
	for(i=0;i<3;++i){
		free(patch_x[i]);
		free(patch_z[i]);
	}
	x = (double *)calloc(N,sizeof(double));
	z = (double *)calloc(N,sizeof(double));
	bin = (unsigned int *)calloc(N,sizeof(unsigned int));
	for(i=0;i<3;++i){
		patch_x[i]=(double *)calloc(N,sizeof(double));
		patch_z[i]=(double *)calloc(N,sizeof(double));
	}
	patches = (unsigned short int*)calloc(N,sizeof(unsigned short int));
	current_particle = 0;
	return SKIP_ONE;
}

file_state read_particle(const char *zeile){
	char species_indicator;
	double a;
	int vars=0;
	if((vars=sscanf(zeile,"%[CN]\t%lf\t%lf\t%lf",&species_indicator,&(x[current_particle]),&(z[current_particle]),&a)) != 4){
		printf("Read %d vars instead of 4 (species,x,z,a).\n%s\n",vars,zeile);
		return ERROR;
	}
	bin[current_particle]=binof(z[current_particle]);
	switch(species_indicator){
		case 'N':
			patches[current_particle]=2;
			patch_x[0][current_particle]=x[current_particle]+sigma/2.0*cos(a);
			patch_z[0][current_particle]=z[current_particle]+sigma/2.0*sin(a);
			patch_x[1][current_particle]=x[current_particle]+sigma/2.0*cos(a+M_PI);
			patch_z[1][current_particle]=z[current_particle]+sigma/2.0*sin(a+M_PI);
			break;
		case 'C':
			patches[current_particle]=3;
			patch_x[0][current_particle]=x[current_particle]+sigma/2.0*cos(a);
			patch_z[0][current_particle]=z[current_particle]+sigma/2.0*sin(a);
			patch_x[1][current_particle]=x[current_particle]+sigma/2.0*cos(a+2.0/3.0*M_PI);
			patch_z[1][current_particle]=z[current_particle]+sigma/2.0*sin(a+2.0/3.0*M_PI);
			patch_x[2][current_particle]=x[current_particle]+sigma/2.0*cos(a+4.0/3.0*M_PI);
			patch_z[2][current_particle]=z[current_particle]+sigma/2.0*sin(a+4.0/3.0*M_PI);
			break;
		default:
			printf("Error. Unknown species.\n");
			return ERROR;
	}
	bin_patches[bin[current_particle]]+=patches[current_particle];
	bin_bonds[bin[current_particle]]+=bondsof();
	bin_n[bin[current_particle]]+=1;
	if(patches[current_particle]==3){
		bin_composition[bin[current_particle]]+=1;
	}
	if(++current_particle==N) return END_OF_FRAME;
	else return READ_PARTICLE;
}

unsigned int binof(double z){
	if(z>HEIGHT){
		printf("Error. Invalid point.\n");
		exit(1);
	}
	return (unsigned int)(floor(z/BINHEIGHT));
}

unsigned int bondsof(void){
	unsigned int bonds=0;
	int j=0;
	for(j=0;j<current_particle;++j){
		if(hasbond(current_particle,j)){
			bin_bonds[bin[j]]+=1;
			++bonds;
		}
	}
	return bonds;
}

bool hasbond(int p1, int p2){
	if(distance(p1,p2)>sigma+delta) return false;
	int i=0,j=0;
	for(i=0;i<patches[p1];++i){
		for(j=0;j<patches[p2];++j){
			if(patch_distance(p1,i,p2,j)<=delta) return true;
		}
	}
	return false;
}

double patch_distance(int p1, int patch1, int p2, int patch2){
	return sqrt(pow(patch_x[patch1][p1]-patch_x[patch2][p2],2)+pow(patch_z[patch1][p1]-patch_z[patch2][p2],2));
}

double distance(int p1, int p2){
	return sqrt(pow(x[p1]-x[p2],2)+pow(z[p1]-z[p2],2));
}
