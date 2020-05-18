#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include "mycustom/objective.h"

void readparameter(int* N, int* D, int* UPBD, int* LWBD, int* MAXITER, int* k, double* WMIN, double* WMAX, double* C1, double* C2) {
	FILE* filePtr;
	char filename[] = "rmo.txt";
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
			fscanf(filePtr, "%s%d", &str, &*k);
			fscanf(filePtr, "%s%f", &str, &*WMIN);
			fscanf(filePtr, "%s%f", &str, &*WMAX);
			fscanf(filePtr, "%s%f", &str, &*C1);
			fscanf(filePtr, "%s%f", &str, &*C2);
		}
	}
}

double randgen() {
	return rand()/(double)RAND_MAX;
}

void INIT(int N, int D, double PAR[][D], int UPBD, int LWBD) {
	int i, j;
	for(i=0; i<N; i++) {
		for(j=0; j<D; j++) {
			PAR[i][j] = randgen()*(UPBD-LWBD) + LWBD;
		}
	}
}

void VELOCITY(int N, int D, double VMAX, double VEL[][D]) {
	int i, j;
	for(i=0; i<N; i++) {
		for(j=0; j<D; j++) {
//			VEL[i][j] = randgen()*VMAX;
			VEL[i][j] = (randgen()*2-1)*VMAX;
		}
	}
}

void SPRINKLE(int N, int D, double WEIGHT, double VEL[][D], double PAR[][D], double CENTER[]) {
	int i, j;
	for(i=0; i<N; i++) {
		for(j=0; j<D; j++) {
			PAR[i][j] = VEL[i][j]*WEIGHT + CENTER[j];
		}
	}
}

void SPRINKLE2(int N, int D, double WEIGHT, double VMAX, double PAR[][D], double CENTER[]) {
	int i, j;
	for(i=0; i<N; i++) {
		for(j=0; j<D; j++) {
			PAR[i][j] = (randgen()*2-1)*VMAX*WEIGHT + CENTER[j];
		}
	}
}

void SETBEST(int N, int D, double PAR[][D], double FIT[], double CENTER[], double GBESTPOS[], double* GBESTFIT) {
	int i, j, minId=0;
	for(i=1; i<N; i++) {
		if(FIT[i] < FIT[minId]) {
			minId = i;
		}
	}
	for(j=0; j<D; j++) {
		GBESTPOS[j] = PAR[minId][j];
		CENTER[j] = PAR[minId][j];
	}
	*GBESTFIT = FIT[minId];
}

void RBEST(int N, int D, double PAR[][D], double FIT[], double RBESTPOS[], double* RBESTFIT) {
	int i, j, minId=0;
	double BESTFIT = FIT[minId];
	for(i=1; i<N; i++) {
		if(FIT[i] < BESTFIT) {
			minId=i;
			BESTFIT = FIT[minId];
		}
	}
	for(j=0; j<D; j++) {
		RBESTPOS[j] = PAR[minId][j];
	}
	*RBESTFIT = FIT[minId];	
}

void GBEST(int D, double RBESTPOS[], double RBESTFIT, double GBESTPOS[], double* GBESTFIT) {
	if(RBESTFIT < *GBESTFIT) {
		*GBESTFIT = RBESTFIT;
		int j;
		for(j=0; j<D; j++) {
			GBESTPOS[j] = RBESTPOS[j];
		}
	}
}

void THECENTER(int D, double C1, double C2, double CENTER[], double RBESTPOS[], double GBESTPOS[]) {
	int j;
	double UPDATE=0;
	for(j=0; j<D; j++) {
		UPDATE = C1*(GBESTPOS[j] - CENTER[j]) + C2*(RBESTPOS[j] - CENTER[j]);
		CENTER[j] = CENTER[j] + UPDATE;
	}
}

void AMMEND(int UPBD, int LWBD, int N, int D, double PAR[][D]) {
	int i, j;
	for(i=0; i<N; i++) {
		for(j=0; j<D; j++) {
			if(PAR[i][j] > (double)UPBD) {
				PAR[i][j] = (double)UPBD;
			}
			if(PAR[i][j] < (double)LWBD) {
				PAR[i][j] = (double)LWBD;
			}
		}
	}
}

void printAll(int N, int D, double P[][D]) {
	int i, j;
	for(i=0; i<N; i++) {
		printf("P%d:\t", i+1);
		for(j=0; j<D; j++) {
			printf("%f  ",P[i][j]);
		}
		printf("\n");
	}
}

int main(int argc, const char * argv[]) {
	int N=50, D=10, UPBD=100, LWBD=-100, MAXITER=500, k=100;
//	double WEIGHT, WMIN=0.0, WMAX=1.0, VMAX=(UPBD-LWBD)/k, C1=0.8, C2=0.9;
	double WEIGHT, WMIN=0.0, WMAX=1.0, VMAX=1, C1=0.8, C2=0.9;
	readparameter(&N,&D,&UPBD,&LWBD,&MAXITER,&k,&WMIN,&WMAX,&C1,&C2);
//	printf("N=%d D=%d UP=%d LW=%d T=%d k=%d %f %f %f %f\n", N, D, UPBD, LWBD, MAXITER, k, WMIN, WMAX, C1, C2);

	VMAX = (double)(UPBD-LWBD)/k;

	double PAR[N][D];
	double VEL[N][D];
	
	double FIT[N];
	
	double RBESTPOS[D];
	double RBESTFIT;
	
	double GBESTPOS[D];
	double GBESTFIT;
	
	double CENTER[D];
	
	srand( time(NULL) );
	
	INIT(N,D,PAR,UPBD,LWBD);
	F1(N,D,PAR,FIT);
	SETBEST(N,D,PAR,FIT,CENTER,GBESTPOS,&GBESTFIT);
//	printAll(N,D,PAR);
	printf("RADIAL MOVEMENT OPTIMIZATION (RMO)\nRUNNING...\nINITIAL BEST: %f\n", GBESTFIT);
	int t;
	for(t=0; t<MAXITER; t++) {
//		VELOCITY(N,D,VMAX,VEL);
//		printAll(N,D,VEL);
		WEIGHT = WMAX - (WMAX-WMIN)*(t/(double)MAXITER);
//		SPRINKLE(N,D,WEIGHT,VEL,PAR,CENTER);
		SPRINKLE2(N,D,WEIGHT,VMAX,PAR,CENTER);
		AMMEND(UPBD,LWBD,N,D,PAR);
		F1(N,D,PAR,FIT);
		RBEST(N,D,PAR,FIT,RBESTPOS,&RBESTFIT);
		THECENTER(D,C1,C2,CENTER,RBESTPOS,GBESTPOS);
		GBEST(D,RBESTPOS,RBESTFIT,GBESTPOS,&GBESTFIT);
//		printf("%f - %f\n",RBESTFIT,GBESTFIT);
//		printf("%f\n",WEIGHT);
	}
	printf("FINAL BEST: %f\n", GBESTFIT);
	printf("DONE...\n=====================================\n");
	getch();
	return 0;
}

