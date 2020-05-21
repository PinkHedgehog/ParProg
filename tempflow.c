#include <stdio.h>
#include <stdlib.h>

#define L 1.0
#define T1 1.0
#define T2 2.0

int main(int argc, char **argv)
{
	if (argc != 3) {
		printf("Usage: %s T M\n", argv[0]);
		return 1;
	}
	//double start = MPI_Wtime();
	double T = atof(argv[1]);
	int M = atoi(argv[2]);
	double h = L / M;
	double tau = 0.3 * h * h;
	int N = T / tau;
	double *u0 = (double*) malloc(sizeof(double) * M);
	double *u1 = (double*) malloc(sizeof(double) * M);
	int m, n;
	
	/* Задаем начальные условия. */
	for (m = 0; m < M; m++) {
		u0[m] = u1[m] = 0.0;
	}
	/* Задаем граничные условия. */
	u0[0] = u1[0] = T1;
	u0[M - 1] = u1[M - 1] = T2;
	
	/* Интегрируем по времени. */
	for (n = 0; n < N; n++) {
		for (m = 1; m < M - 1; m++) {
			u1[m] = u0[m] + tau / h / h  * (u0[m-1] - 2.0 * u0[m] + u0[m+1]);
		}
		double *t = u0;
		u0 = u1;
		u1 = t;
	}
	
	/* Вывод на экран. */
	for (m = 0; m < M; m++) {
		printf("%lf %lf\n", m * h, u1[m]);
	}
	double end = MPI_Wtime();
	printf("%lf\n", end - start);

	free(u0);
	free(u1);
	return 0;
}
