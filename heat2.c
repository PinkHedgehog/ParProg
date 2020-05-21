#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define Length 1.0
#define Temperature_1 1.0
#define Temperature_2 5.0

int main(int argc, char **argv)
{
    if (argc < 3)
	{
        printf("Usage: %s <time> <step> \n", argv[0]);
		exit(1);
	}


	int size, rank;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	// Target time
	double Time = atof(argv[1]);
	if (Time < 0) {
		if(rank == 0){
			printf("Sorry, timemachine hasn't been invented yet!");
			MPI_Finalize();
			return EXIT_FAILURE;
		}
	}

	int M = atoi(argv[2]);
	if (M < 2) {
		if(rank == 0){
			printf("Invalid values!\n");
			MPI_Finalize();
			return EXIT_FAILURE;
		}
	} else if(M <= size) {
		if(rank == 0){
			printf("Incorrect number of partitions\n");
			MPI_Finalize();
			return EXIT_FAILURE;
		}
	}

	double h = Length / M;
	double tau = 0.3 * h * h;
	int N = Time / tau;

    // T(n), T(n+1)
	double *u0 = (double*) malloc(sizeof(double) * M);
	double *u1 = (double*) malloc(sizeof(double) * M);

	size_t m, n;

	double time = MPI_Wtime();
	
	// At the beginning (f(x) = 0 )
	for (m = 0; m < M; m++) 
    {
		u0[m] = u1[m] = 0.0;
	}

	// Edges
	u0[0] = u1[0] = Temperature_1;
	u0[M - 1] = u1[M - 1] = Temperature_2;
	
	size_t *left_index = (size_t*) malloc(sizeof(size_t) * size + 1);
	left_index[0] = 1;
	left_index[size] = M - 1;
	for(int i = 1; i < size; i++) {
		left_index[i] = left_index[i - 1] + (M / size) + ((i - 1) < ((M % size) - 2));
	}

	 // Цикл по времени
	for (n = 0; n < N; n++) {
		// Swap
		if(rank != 0) {
			MPI_Send(u0 + left_index[rank], 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
			MPI_Recv(u0 + left_index[rank] - 1, 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
		if(rank != size - 1) {
			MPI_Send(u0 + left_index[rank + 1] - 1, 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
			MPI_Recv(u0 + left_index[rank + 1], 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}

		// "Main"
		for (m = left_index[rank]; m < left_index[rank + 1]; m++) {
			u1[m] = u0[m] + 0.3  * (u0[m - 1] - 2.0 * u0[m] + u0[m + 1]);
		}
		// Refresh
		double *t = u0;
		u0 = u1;
		u1 = t;
	}

	if(size > 1) {
		if(rank == 0) {
			for(int i = 1; i < size; i++) {
				MPI_Recv(u1 + left_index[i], left_index[i + 1] - left_index[i], MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			}
		} else {
			MPI_Send(u1 + left_index[rank], left_index[rank + 1] - left_index[rank], MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
		}
	}

	// Рассчитываем время работы программы
	time = MPI_Wtime() - time;
	
	if(rank == 0) {
		/*for (m = 0; m < M; m++) {
			printf("%lf %lf\n", m * h, u1[m]);
		}*/
		printf("%d %lf\n", size, time);
	}
	
	free(u0);
	free(u1);

    MPI_Finalize();
	return EXIT_SUCCESS;
}
