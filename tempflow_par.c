#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define L 1.0
#define T1 1.0
#define T2 2.0

int main(int argc, char **argv)
{
	if (argc != 3) {
		printf("Usage: %s T M\n", argv[0]);
		return 1;
	}
	
	MPI_Init(&argc, &argv)
	double start = MPI_Wtime();
	double T = atof(argv[1]);
	int M = atoi(argv[2]);
	double h = L / M;
	double tau = 0.3 * h * h;
	int N = T / tau;
	double *u0 = (double*) malloc(sizeof(double) * M);
	double *u1 = (double*) malloc(sizeof(double) * M);
	int m, n;
	
	/* начало и конец зоны ответсвенности текущего процесса. */
	int sk, fk;
	
	/* Число процессов. */
	int p;
	/* Номер текущего процесса (ранк). */
	int k;
	MPI_Comm_rank(MPI_COMM_WORLD, &k);
	MPI_Comm_size(MPI_COMM_WORLD, &p);
	
	int mk = M / p;
	/*
	sk = k * mk + ((k < M % p) ? k : M % p);
	fk = sk + mk + (k < M % p);
	*/
	sk = k * mk;
	fk = sk + mk;
	if (k == 0) sk = 1;
	if (k == p - 1) fk = M - 1;
	
	//printf("#%d: sk = %d, fk = %d\n", k, sk, fk);
	
	/* Задаем начальные условия. */
	for (m = 0; m < M; m++) {
		u0[m] = u1[m] = 0.0;
	}
	/* Задаем граничные условия. */
	u0[0] = u1[0] = T1;
	u0[M - 1] = u1[M - 1] = T2;
	
	/* Интегрируем по времени. */
	for (n = 0; n < N; n++) {
		for (m = sk; m < fk; m++) {
			u1[m] = u0[m] + tau / h / h  * (u0[m-1] - 2.0 * u0[m] + u0[m+1]);
			//u1[m] = k;
		}
		double *t = u0;
		u0 = u1;
		u1 = t;
	}
	
	/* Сбор данных. */
	if (k == 0) {
		int i;
		for (i = 1; i < p; i++) {
			int mi = mk;
			if (i == p - 1) mi = M - 1 - i * mk;
			MPI_Recv(u1 + i * mk, mi, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
	} else {
		MPI_Send(u1 + sk, fk - sk, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	}
	/* Вывод на экран. */
	// if (k == 0) {
	// 	for (m = 0; m < M; m++) {
	// 		printf("%lf %lf\n", m * h, u1[m]);
	// 	}
	// }
	
	free(u0);
	free(u1);
	if (k == 0)
	{
		double end = MPI_Wtime();
		printf("%lf\n", end - start);
	}
	MPI_Finalize();
	

	return 0;
}
