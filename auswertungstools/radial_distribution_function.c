#include <stdio.h>
#include <math.h>
#include <stdlib.h>

typedef enum {READNUM,SKIPLINE,READPOS,END} state;

const double binwidth = 0.2;
const double height = 100;
const double width = 100;
double *density;
double cutoff = 0;
int N = 0;
state file_state = READNUM;
double *x=NULL;
double *z=NULL;
int frames = 0;

double min(double a, double b){
	return (a<b?a:b);
}

int bin_index(double r){
	return (int)(floor(r/binwidth));
}

int main(int argc, char **argv){
	if(argc > 1){
		cutoff = min(height,width)/2.0;
		int bins=(int)ceil(cutoff/binwidth);
		density = (double *)calloc(bins,sizeof(double));
		int i=0,j=0;
		FILE *datei=fopen(argv[1],"r");
		size_t len = 50;
		char *zeile=calloc(50,sizeof(char));
		char species;
		double discard;
		double dx,dz,r;
		while(getline(&zeile,&len,datei) != -1){
			switch(file_state){
				case READNUM:
					sscanf(zeile,"%d",&N);
					i=0;
					free(x);
					free(z);
					x = (double*)calloc(N,sizeof(double));
					z = (double*)calloc(N,sizeof(double));
					++frames;
					file_state=SKIPLINE;
					break;
				case SKIPLINE:
					file_state = READPOS;
					break;
				case READPOS:
					sscanf(zeile,"%[CN] %lf %lf %lf",&species,&x[i],&z[i],&discard);
					for(j=0;j<i;++j){
						dx=x[i]-x[j];
						if(abs(dx)<=cutoff){
							dz=z[i]-z[j];
							if(abs(dz)<=cutoff){
								r=sqrt(dx*dx+dz*dz);
								if(r <= cutoff){
									density[bin_index(r)]++;
								}
							}
						}
					}
					++i;
					if(i==N){
						printf("read frame %d.\n",frames);
						file_state=READNUM;
					}
					break;
				default:
					printf("An unexpected thing occured\n");
					break;
			}
		}
		free(zeile);
		fclose(datei);
		double norm;
		for(i=0;i<bins;++i){
			norm=M_PI*(2.0*i+1.0)*binwidth*binwidth*((double)N)/(height*width);
			density[i]/=norm*frames;
		}

		FILE *output=fopen("rdf.dat","w");
		fprintf(output,"#bin\tr}trdf(r)\n");
		for(i=0;i<bins;++i){
			fprintf(output,"%d\t%.2f\t%.5e\n",i,(i+0.5)*binwidth,density[i]);
		}
		free(density);
		fclose(output);
	}else{
		printf("Give me a filename!\n");
	}
}
