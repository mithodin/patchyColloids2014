#include <stdio.h>
#include <math.h>
#include <stdlib.h>

typedef enum {READNUM,SKIPLINE,READPOS,END} state;
typedef enum {THREEPATCH,TWOPATCH} species;

const double binheight = 1.0;
const double height = 200;
const double width = 100;
double *pbond;
int N = 0;
state file_state = READNUM;
double *x=NULL;
double *z=NULL;
double *a=NULL;
species *sp=NULL;
int *totalpatches = NULL;
int frames = 0;
const double sigma = 1.0;
const double delta = 0.11965683746373795115;

int bin_index(double z){
	return (int)(floor(z/binheight));
}

int patches(int i){
	switch(sp[i]){
		case THREEPATCH:
			return 3;
			break;
		case TWOPATCH:
			return 2;
			break;
		default:
			printf("Error!\n");
			exit(1);
	}
}

double patchPositionX(int k, int i){
	if( sp[k] == THREEPATCH ){
		return x[k] + sigma/2.0*cos(a[k] + (i%3)*2.0/3.0*M_PI);
	}else if( sp[k] == TWOPATCH ){
		return x[k] + sigma/2.0*cos(a[k] + (i%2)*M_PI);
	}else{
		printf("invalid patch position.\n");
		exit(-1);
	}
}

double patchPositionZ(int k, int i){
	if( sp[k] == THREEPATCH ){
		return z[k] + sigma/2.0*sin(a[k] + (i%3)*2.0/3.0*M_PI);
	}else if( sp[k] == TWOPATCH ){
		return z[k] + sigma/2.0*sin(a[k] + (i%2)*M_PI);
	}else{
		printf("invalid patch position.\n");
		exit(-1);
	}
}

int bonds(int k, int l){
	int i = 0;
	int j = 0;
	double d = sqrt(pow(x[k]-x[l],2)+pow(z[k]-z[l],2));
	if(d > sigma+delta) return 0;
	for(i = 0; i < patches(k); i++){
		for( j = 0; j < patches(l); ++j ){
			double p1X = patchPositionX(k,i);
			double p1Z = patchPositionZ(k,i);
			double p2X = patchPositionX(l,j);
			double p2Z = patchPositionZ(l,j);
				
			d = sqrt((p1X-p2X)*(p1X-p2X)+(p1Z-p2Z)*(p1Z-p2Z));
			if( d <= delta ){
				return 1;
			}
		}
	}
	return 0;
}

int main(int argc, char **argv){
	if(argc > 1){
		int bins=(int)ceil(height/binheight);
		pbond = (double *)calloc(bins,sizeof(double));
		totalpatches = (int *)calloc(bins,sizeof(int));
		int i=0,j=0;
		FILE *datei=fopen(argv[1],"r");
		size_t len = 80;
		char *zeile=calloc(80,sizeof(char));
		char sp0;
		double dx,dz;
		while(getline(&zeile,&len,datei) != -1){
			switch(file_state){
				case READNUM:
					sscanf(zeile,"%d",&N);
					i=0;
					free(x);
					free(z);
					free(a);
					free(sp);
					x = (double*)calloc(N,sizeof(double));
					z = (double*)calloc(N,sizeof(double));
					a = (double*)calloc(N,sizeof(double));
					sp = (species*)calloc(N,sizeof(species));
					file_state=SKIPLINE;
					break;
				case SKIPLINE:
					file_state = READPOS;
					break;
				case READPOS:
					sscanf(zeile,"%[CN] %lf %lf %lf",&sp0,&x[i],&z[i],&a[i]);
					switch(sp0){
						case 'C':
							sp[i]=THREEPATCH;
							break;
						case 'N':
							sp[i]=TWOPATCH;
							break;
						default:
							printf("Error! Wrong Species!\n");
							exit(1);
					}
					totalpatches[bin_index(z[i])]+=patches(sp[i]);
					for(j=0;j<i;++j){
						pbond[bin_index(z[i])]+=bonds(i,j);
					}
					++i;
					if(i==N){
						frames++;
						printf("read frame %d.\n",frames);
						file_state=READNUM;
					}
					break;
				default:
					printf("An unexpected thing occured\n");
					break;
			}
		}
		//free(zeile);
		fclose(datei);
		for(i=0;i<bins;++i){
			pbond[i]/=totalpatches[i];
		}

		FILE *output=fopen("pbond.dat","w");
		fprintf(output,"#bin\tz\tpbond(z)\n");
		for(i=0;i<bins;++i){
			fprintf(output,"%d\t%.2f\t%.5e\n",i,(i+0.5)*binheight,pbond[i]);
		}
		fclose(output);
		//free(pbond);
		free(totalpatches);
	}else{
		printf("Give me a filename!\n");
	}
}
