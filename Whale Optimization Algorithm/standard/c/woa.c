#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include "mycustom/objective.h"

#define PI 3.1415926535897931

void readparameter(int* N, int* D, int* UPBD, int* LWBD, int* MAXITER, double* b) {
	FILE* filePtr;
	char filename[] = "woa.txt";
	char str[50];
	if((filePtr = fopen(filename,"r")) == NULL) {
		printf("FILE COULD NOT BE OPENED.\n");
		exit(0);
	}
	else {
		while(!feof(filePtr)) {
			fscanf(filePtr, "%s%d", &str, &*N); //printf("%d\n",*N);
			fscanf(filePtr, "%s%d", &str, &*D);
			fscanf(filePtr, "%s%d", &str, &*MAXITER);
			fscanf(filePtr, "%s%d", &str, &*UPBD);
			fscanf(filePtr, "%s%d", &str, &*LWBD);
			fscanf(filePtr, "%s%f", &str, &*b);
		}
	}
}

double randgen() {
	return rand()/(double)RAND_MAX;
}

void INIT(int N, int D, double W[][D], int UPBD, int LWBD) {
	int i, j;
	for(i=0; i<N; i++) {
		for(j=0; j<D; j++) {
			W[i][j] = randgen()*(UPBD-LWBD) + LWBD;
		}
	}
}

void SETBEST(int N, int D, double W[][D], double FIT[], double GBESTPOS[], double* GBESTFIT) {
	int i, j, minId=0;
	*GBESTFIT = FIT[minId];
	for(i=1; i<N; i++) {
		if(FIT[i] < *GBESTFIT) {
			minId=i;
			*GBESTFIT = FIT[minId];
		}
	}
	for(j=0; j<D; j++) {
		GBESTPOS[j] = W[minId][j];
	}
}

void SHRINKINGCIRCLING(int D, double AA, double CC, double W[], double GBESTPOS[]) {
	int j;
	for(j=0; j<D; j++) {
		W[j] = GBESTPOS[j] - AA*(fabs(CC*GBESTPOS[j] - W[j]));
	}
}

void SEARCHFORPREY(int D, double AA, double CC, double W[], double RANDWHALE[]) {
	int j;
	double distance;
	for(j=0; j<D; j++) {
		distance = fabs(CC*RANDWHALE[j] - W[j]);
		W[j] = RANDWHALE[j] - AA*distance;
	}
}

void SPIRALUPDATING(int D, double AA, double CC, double b, double l, double W[], double GBESTPOS[]) {
	int j;
	double distance;
	for(j=0; j<D; j++) {
		distance = fabs(GBESTPOS[j] - W[j]);
		W[j] = distance * exp(b*l) * cos(2*PI*l) + GBESTPOS[j];
	}
}

void GBEST(int N, int D, double W[][D], double FIT[], double GBESTPOS[], double *GBESTFIT) {
	int i, j, minId=0;
	*GBESTFIT = FIT[minId];
	for(i=1; i<N; i++) {
		if(FIT[i] < *GBESTFIT) {
			minId = i;
			*GBESTFIT = FIT[minId];
		}
	}
	for(j=0; j<D; j++) {
		GBESTPOS[j] = W[minId][j];
	}
}

void AMMEND(int UPBD, int LWBD, int N, int D, double W[][D]) {
	int i, j;
	for(i=0; i<N; i++) {
		for(j=0; j<D; j++) {
			if(W[i][j] > (double)UPBD) {
				W[i][j] = (double)UPBD;
			}
			if(W[i][j] < (double)LWBD) {
				W[i][j] = (double)LWBD;
			}
		}
	}
}

int main(int argc, const char * argv[]) {
	
	int N=30, D=10, MAXITER=100, UPBD=100, LWBD=-100;
	double b=1.5;
	readparameter(&N,&D,&UPBD,&LWBD,&MAXITER,&b);
	srand(time(NULL));
	double WHALE[N][D];
	double FIT[N];
	double GBESTPOS[D];
	double GBESTFIT;
	
	INIT(N,D,WHALE,UPBD,LWBD);
	F1(N,D,WHALE,FIT);
	SETBEST(N,D,WHALE,FIT,GBESTPOS,&GBESTFIT);
	printf("WHALE OPTIMIZATION ALGORITHM (WOA)\nRUNNING...\nINITIAL BEST: %f\n", GBESTFIT);
	int t, w, randId;
	double a, a2, l, AA, CC, prob;
	for(t=0; t<MAXITER; t++) {
		//a decreases linearly fron 2 to 0 in Eq. (2.3)
    	a=2 - t*((2-0)/MAXITER);
		for(w=0; w<N; w++) {
			AA = randgen()*(2*a)-a;
    		CC = randgen()*2;
    		prob = randgen();
    		if(prob < 0.5) {
    			if(fabs(AA) < 1) {
    				SHRINKINGCIRCLING(D,AA,CC,WHALE[w],GBESTPOS);
				}
				else {	//fabs(AA) >= 1
					randId = rand()%N; // index of a random whale
					SEARCHFORPREY(D,AA,CC,WHALE[w],WHALE[randId]);
				}
			}
			else {	// prob >= 0.5
				l = randgen()*2-1;
				//a2=-1 - t*((-1-(-2))/MAXITER);
				//l = (a2-1)*randgen() + 1;
				SPIRALUPDATING(D,AA,CC,b,l,WHALE[w],GBESTPOS);
			}
		}
		AMMEND(UPBD,LWBD,N,D,WHALE);
		F1(N,D,WHALE,FIT);
		GBEST(N,D,WHALE,FIT,GBESTPOS,&GBESTFIT);
//		printf("ITER %d %f\n",t, GBESTFIT);
	}
	
	printf("FINAL BEST: %f\n", GBESTFIT);
	printf("DONE...\n=====================================\n");
//	
//	int m;
//	for(m=0; m<D; m++) {
//		printf("%f\n", GBESTPOS[m]);
//	}
	
	getch();
	return 0;
}

