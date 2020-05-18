#include <stdio.h>
#include <math.h>

void F1(int N, int D, double P[][D], double F[]) {
	int i, j;
	double fitness;
	for(i=0; i<N; i++) {
		fitness=0;
		for(j=0; j<D; j++) {
			fitness += pow(P[i][j],2);
		}
		F[i] = fitness;
	}
}

