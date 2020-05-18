#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include "mycustom/objective.h"

void readparameter(int* N, int* D, int* UPBD, int* LWBD, int* MAXITER, double* WEIGHT, double* C1, double* C2) {
	FILE* filePtr;
	char filename[] = "pso.txt";
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
			fscanf(filePtr, "%s%f", &str, &*WEIGHT);
			fscanf(filePtr, "%s%f", &str, &*C1);
			fscanf(filePtr, "%s%f", &str, &*C2);
		}
	}
}

double randgen() {
	return rand()/(double)RAND_MAX;
}

void INIT(int N, int D, double SWARM[][D], int UPBD, int LWBD) {
	int i, j;
	for(i=0; i<N; i++) {
		for(j=0; j<D; j++) {
			SWARM[i][j] = randgen()*(UPBD-LWBD) + LWBD;
		}
	}
}

void SETPBEST(int N, int D, double SWARM[][D], double PBESTPOS[][D], double FIT[], double PBESTFIT[]) {
	int i, j;
	for(i=0; i<N; i++) {
		for(j=0; j<D; j++) {
			PBESTPOS[i][j] = SWARM[i][j];
		}
		PBESTFIT[i] = FIT[i];
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

void printFit(int N, double FIT[], double PBESTFIT[]) {
//	printf("\n");
	int i;
	for(i=0; i<N; i++) {
		printf("P[%d] %f   %f\n", i, FIT[i], PBESTFIT[i]);
	}
}

void SETBEST(int N, int D, double SWARM[][D], double FIT[], double GBESTPOS[], double* GBESTFIT) {
	int i, j, minId=0;
	for(i=1; i<N; i++) {
		if(FIT[i] < FIT[minId]) {
			minId = i;
		}
	}
	for(j=0; j<D; j++) {
		GBESTPOS[j] = SWARM[minId][j];
	}
	*GBESTFIT = FIT[minId];
}

void GBEST(int N, int D, double PBESTPOS[][D], double PBESTFIT[], double GBESTPOS[], double* GBESTFIT) {
	int i, j;
	for(i=0; i<N; i++) {
		if(PBESTFIT[i] < *GBESTFIT) {
			*GBESTFIT = PBESTFIT[i];
			for(j=0; j<D; j++) {
				GBESTPOS[j] = PBESTPOS[i][j];
			}
		}
	}
}

void GETPBEST(int N, int D, double SWARM[][D], double PBESTPOS[][D], double FIT[], double PBESTFIT[]) {
	int i, j;
	for(i=0; i<N; i++) {
		if(FIT[i] < PBESTFIT[i]) {
			PBESTFIT[i] = FIT[i];
			for(j=0; j<D; j++) {
				PBESTPOS[i][j] = SWARM[i][j];
			}
		}
	}
}

void VELOCITY(int N, int D, double WEIGHT, double C1, double C2, double VEL[][D], double SWARM[][D], double PBESTPOS[][D], double GBESTPOS[]) {
	int i, j;
	for(i=0; i<N; i++) {
		for(j=0; j<D; j++) {
			VEL[i][j] = WEIGHT*VEL[i][j] + C1*(randgen()*(PBESTPOS[i][j] - SWARM[i][j])) + C2*(randgen()*(GBESTPOS[j] - SWARM[i][j]));
		}
	}
}

void POSITION(int N, int D, double SWARM[][D], double VEL[][D]) {
	int i, j;
	for(i=0; i<N; i++) {
		for(j=0; j<D; j++) {
			SWARM[i][j] = SWARM[i][j] + VEL[i][j];
		}
	}
}

void AMMEND(int UPBD, int LWBD, int N, int D, double SWARM[][D]) {
	int i, j;
	for(i=0; i<N; i++) {
		for(j=0; j<D; j++) {
			if(SWARM[i][j] > (double)UPBD) {
				SWARM[i][j] = (double)UPBD;
			}
			if(SWARM[i][j] < (double)LWBD) {
				SWARM[i][j] = (double)LWBD;
			}
		}
	}
}

int main(int argc, const char * argv[]) {
//	int run;
//	for(run=0; run<30; run++) {
//	printf("RUN [%d]\n",run+1);
	
	int N=300, D=100, UPBD=100, LWBD=-100, MAXITER=1000;
	double WEIGHT=1.0, C1=0.8, C2=0.9;
	
	readparameter(&N,&D,&UPBD,&LWBD,&MAXITER,&WEIGHT,&C1,&C2);
//	printf("%d %d %d %d %d %f %f %f\n", N, D, UPBD, LWBD, MAXITER, WEIGHT, C1, C2);
	
	double SWARM[N][D];
	double VEL[N][D];
	double PBESTPOS[N][D];
	
	double FIT[N];
	double PBESTFIT[N];
	
	double GBESTPOS[D];
	double GBESTFIT;
	
	srand( time(NULL) );
	
	INIT(N,D,SWARM,UPBD,LWBD);
	F1(N,D,SWARM,FIT);
	SETPBEST(N,D,SWARM,PBESTPOS,FIT,PBESTFIT);
	SETBEST(N,D,SWARM,FIT,GBESTPOS,&GBESTFIT);

//	printAll(N,D,SWARM); printf("\n");
//	printAll(N,D,PBESTPOS);	
//	printFit(N,FIT,PBESTFIT);
//	printf("\n%f\n",GBESTFIT);
	printf("PARTICLE SWARM OPTIMIZATION (PSO) ALGORITHM\nRUNNING...\nINITIAL BEST: %f\n", GBESTFIT);
	int t;
	for(t=0; t<MAXITER; t++) {
		WEIGHT = 1 - (1-0)*(t/(double)MAXITER);
		VELOCITY(N,D,WEIGHT,C1,C2,VEL,SWARM,PBESTPOS,GBESTPOS);
		POSITION(N,D,SWARM,VEL);
		AMMEND(UPBD,LWBD,N,D,SWARM);
		F1(N,D,SWARM,FIT);
		GETPBEST(N,D,SWARM,PBESTPOS,FIT,PBESTFIT);
		GBEST(N,D,PBESTPOS,PBESTFIT,GBESTPOS,&GBESTFIT);
//		printf("ITER[%d]  %f\n", t+1, GBESTFIT);
//		printf("ITER [%d]\n",t+1);
//		printFit(N,FIT,PBESTFIT);
	}
	printf("FINAL BEST: %f\n", GBESTFIT);
	printf("DONE...\n=====================================\n");
//	}
	getch();
	return 0;
}
