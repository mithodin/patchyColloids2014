#include <stdio.h>
#include <math.h>
#include <stdlib.h>

typedef enum {READNUM,SKIPLINE,READPOS,END} state;
typedef enum {THREEPATCH,TWOPATCH} species;

const double binheight = 2.0;
const double height = 200;
const double width = 100;
double *pbond;
int N = 0;
double *x=NULL;
double *z=NULL;
double *a=NULL;
int *bin=NULL;
species *sp=NULL;
int *totalpatches = NULL;
const double sigma = 1.0;
const double delta = 0.11965683746373795115;

int bin_index(double z){
	return (int)(floor(z/binheight));
}

double distance(double x1, double z1, double x2, double z2){
	return sqrt(pow(x1-x2,2)+pow(z1-z2,2));
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
			printf("Error! sp[i] is not valid!\n");
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
	double d = distance(x[k],z[k],x[l],z[l]);
	if(d > sigma+delta){
		return 0;
	}else if(d < sigma){
		printf("too close! %f\n",d);
	}else{
		for(i = 0; i < patches(k); i++){
			for( j = 0; j < patches(l); ++j ){
				double p1X = patchPositionX(k,i);
				double p1Z = patchPositionZ(k,i);
				double p2X = patchPositionX(l,j);
				double p2Z = patchPositionZ(l,j);
					
				d = distance(p1X,p1Z,p2X,p2Z);
				if( d <= delta ){
					return 1;
				}
			}
		}
	}
	return 0;
}

int protected_index(int max, int index){
	if(index>=max){
		printf("index %d out of bounds!\n",index);
		exit(1);
	}else{
		return index;
	}
}

int main(int argc, char **argv){
	if(argc > 1){
		int bins=(int)ceil(height/binheight);
		pbond = (double *)calloc(bins,sizeof(double));
		totalpatches = (int *)calloc(bins,sizeof(int));
		int i=0,j=0,tmp=0;
		FILE *datei=fopen(argv[1],"r");
		size_t len = 60;
		char *zeile=malloc(60*sizeof(char));
		char sp0;
		state file_state = READNUM;
		int frames = 0;
		while((sp0=fgetc(datei))!=EOF){
		//while(getline(&zeile,&len,datei) != -1){
			ungetc(sp0,datei);
			switch(file_state){
				case READNUM:
					tmp=fscanf(datei,"%d",&N);
					if(tmp!=1){
						printf("error! Could not read N\n");
						exit(1);
					}
					free(bin);
					free(x);
					free(z);
					free(a);
					free(sp);
					bin = (int*)calloc(N,sizeof(int));
					x = (double*)calloc(N,sizeof(double));
					z = (double*)calloc(N,sizeof(double));
					a = (double*)calloc(N,sizeof(double));
					sp = (species*)calloc(N,sizeof(species));
					i=0;
					file_state=SKIPLINE;
					break;
				case SKIPLINE:
					fgets(zeile,len,datei);
					file_state = READPOS;
					break;
				case READPOS:
					tmp=fscanf(datei,"%[CN] %lf %lf %lf",&sp0,&(x[i]),&(z[i]),&(a[i]));
					if(tmp!=4){
						printf("error! Only %d variables read by fscanf.\n",tmp);
						exit(1);
					}
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
					bin[i]=bin_index(z[i]);
					totalpatches[bin[i]]+=patches(i);
					if(totalpatches[bin[i]]>10000){
						printf("%d in bin %d\n",totalpatches[bin[i]],bin[i]);
					}else if(totalpatches[bin[i]]==0){
						printf("%d in bin %d\n",totalpatches[bin[i]],bin[i]);
					}
					for(j=0;j<i;++j){
						tmp=bonds(i,j);
						pbond[bin[i]]+=tmp;
						pbond[bin[j]]+=tmp;
						/*if(totalpatches[bin[j]]==0){
							printf("j:%d i:%d\n",j,i);
						} */
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
		fclose(datei);
		//free(zeile);
		for(i=0;i<bins;++i){
			if(totalpatches[i]>0){
				pbond[i]/=totalpatches[i];
			}else{
				if(pbond[i]!=0) printf("Something's weird in bin %d.\n%f bonds, %d patches.\n",i,pbond[i],totalpatches[i]);
			}
		}

		FILE *output=fopen("pbond.dat","w");
		fprintf(output,"#bin\tz\tpbond(z)\n");
		for(i=0;i<bins;++i){
			fprintf(output,"%d\t%.2f\t%.5e\n",i,(i+0.5)*binheight,pbond[i]);
		}
		fclose(output);
		free(pbond);
		free(totalpatches);
	}else{
		printf("Give me a filename!\n");
	}
}
