#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include "mycustom/objective.h"

void readparameter(int* N, int* D, int* UPBD, int* LWBD, int* MAXITER) {
	FILE* filePtr;
	char filename[] = "bha.txt";
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
		}
	}
}

double randgen() {
	return rand()/(double)RAND_MAX;
}

void INIT(int N, int D, double STAR[][D], int UPBD, int LWBD) {
	int i, j;
	for(i=0; i<N; i++) {
		for(j=0; j<D; j++) {
			STAR[i][j] = randgen()*(UPBD-LWBD) + LWBD;
		}
	}
}

void SETBEST(int N, int D, double STAR[][D], double BH[], double GBESTPOS[], double FIT[], double* BHFIT, double* GBESTFIT) {
	int i, j, minId=0;
	*GBESTFIT = FIT[minId];
	for(i=1; i<N; i++) {
		if(FIT[i] < *GBESTFIT) {
			minId = i;
			*GBESTFIT = FIT[minId];
			*BHFIT = FIT[minId];
		}
	}
	
	for(j=0; j<D; j++) {
		GBESTPOS[j] = STAR[minId][j];
		BH[j] = STAR[minId][j];
	}
}

void MOVESTAR(int N, int D, double STAR[][D], double BLACKHOLE[]) {
	int i, j;
	for(i=0; i<N; i++) {
		for(j=0; j<D; j++) {
			STAR[i][j] = STAR[i][j] + randgen()*(BLACKHOLE[j] - STAR[i][j]);
		}
	}
}

void SWAP(int N, int D, double STAR[][D], double* FIT, double BLACKHOLE[], double* BHFIT) {
	int i, j, minId=0;
	double BESTFIT = FIT[minId];
	double TEMP[D];
	double TEMPFIT;
	for(i=1; i<N; i++) {
		if(FIT[i] < BESTFIT) {
			minId = i;
			BESTFIT = FIT[minId];
		}
	}
	if(BESTFIT < *BHFIT) {
		TEMPFIT = *BHFIT;
		*BHFIT = FIT[minId];
		FIT[minId] = *BHFIT;
		
		for(j=0; j<D; j++) {
			TEMP[j] = BLACKHOLE[j];
			BLACKHOLE[j] = STAR[minId][j];
			STAR[minId][j] = TEMP[j];
		}
	}
}

void CROSSINGHORIZON(int N, int D, int UPBD, int LWBD, double R, double STAR[][D], double BH[]) {
	double DISTANCE;
	int i, j, k;
	for(i=0; i<N; i++) {
		DISTANCE = 0;
		for(j=0; j<D; j++) {
			DISTANCE +=  pow(BH[j]-STAR[i][j], 2);
		}
		
		if(sqrt(DISTANCE) < R) {
			for(k=0; k<D; k++) {
				STAR[i][k] = LWBD + randgen()*(UPBD-LWBD);
			}
		}
	}
}

double EVENTHORIZON(int N, double BHFIT, double FIT[]) {
	double AVE=0;
	int i;
	for(i=0; i<N; i++) {
		AVE += FIT[i];
	}
	return BHFIT/AVE;
}

void GBEST(int D, double BLACKHOLE[], double BHFIT, double GBESTPOS[], double *GBESTFIT) {
	int j;
	if(BHFIT < *GBESTFIT) {
		*GBESTFIT = BHFIT;
		for(j=0; j<D; j++) {
			GBESTPOS[j] = BLACKHOLE[j];
		}
	}
}

void AMMEND(int UPBD, int LWBD, int N, int D, double STAR[][D]) {
	int i, j;
	for(i=0; i<N; i++) {
		for(j=0; j<D; j++) {
			if(STAR[i][j] > (double)UPBD) {
				STAR[i][j] = (double)UPBD;
			}
			if(STAR[i][j] < (double)LWBD) {
				STAR[i][j] = (double)LWBD;
			}
		}
	}
}

int main(int argc, const char * argv[]) {
	int N=100, D=10, MAXITER=1000, UPBD=100, LWBD=-100;
	readparameter(&N,&D,&UPBD,&LWBD,&MAXITER);
//	printf("%d %d %d %d %d %d %f %f %f %f\n", N, D, UPBD, LWBD, MAXITER, k, WMIN, WMAX, C1, C2);
	double R;
	srand( time(NULL) );
	
	double STAR[N][D];
	double FIT[N];
	double GBESTPOS[D];
	double GBESTFIT;
	double BLACKHOLE[D];
	double BLACKHOLEFIT;
	
	INIT(N,D,STAR,UPBD,LWBD);
	F1(N,D,STAR,FIT);
	SETBEST(N,D,STAR,BLACKHOLE,GBESTPOS,FIT,&BLACKHOLEFIT,&GBESTFIT);
	
	printf("BLACKHOLE (BH) ALGORITHM\nRUNNING...\nINITIAL BEST: %f\n", GBESTFIT);
	int t;
	for(t=0; t<MAXITER; t++) {
		MOVESTAR(N,D,STAR,BLACKHOLE);
		AMMEND(UPBD,LWBD,N,D,STAR);
		F1(N,D,STAR,FIT);
		SWAP(N,D,STAR,FIT,BLACKHOLE,&BLACKHOLEFIT);
		R = EVENTHORIZON(N,BLACKHOLEFIT,FIT);
		CROSSINGHORIZON(N,D,UPBD,LWBD,R,STAR,BLACKHOLE);
		GBEST(D,BLACKHOLE,BLACKHOLEFIT,GBESTPOS,&GBESTFIT);
	}
	printf("FINAL BEST: %f\n", GBESTFIT);
	printf("DONE...\n=====================================\n");
	getch();
	return 0;
}

