#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <omp.h>

#define Length 1.0
#define Temperature_1 1.0
#define Temperature_2 1.0

int main(int argc, char **argv)
{


	if (argc < 4) // file, Time, Num_Proc, Num_div // ./a.out 1 [1,2,3,4] 100000
	{
        printf("\n\n\nSo few arguments!!! You must get 3 arguments: Time, Num_Proc, Num_div.\n\n\n");
		exit(0);
	}

	double Time = atof(argv[1]); // Считываем момент времени, в который хотим узнать распределение температуры
	if (Time < 0) 
	{
		printf("It is impossible!");
		exit(1);
	}

	int size = atoi(argv[2]);
	if (size <= 0)
	{
		printf("It is impossible");
		exit(1);
	}

	int M = atoi(argv[3]); // Считываем число разбиений по координате
	if (M < 2) 
	{
		printf("Bad values!\n");
		exit(2);
	} 
	else if(M <= size) 
	{ // Плохая мелкость разбиения => не все процессы задействованы
		printf("So many processes to solve this task!\n");
		exit(3);
	}


	double h = Length / M;
	double tau = 0.3 * h * h; // Шаг по времени
	int N = Time / tau;// Число разбиений по времени

    // Массивы температуры для момента времени n и n + 1 соответственно
	double* u0 = (double*) malloc(sizeof(double) * M);
	double* u1 = (double*) malloc(sizeof(double) * M);
	//-----------------------------------------------------------------

	size_t m, n;
	
	// ---Начальные условия ( u(0,x) = 0 )---
	for (m = 0; m < M; m++) 
    {
		u0[m] = u1[m] = 0.0;
	}
	// --------------------------------------

	// --------Граничные условия-------------
	u0[0] = u1[0] = Temperature_1;
	u0[M - 1] = u1[M - 1] = Temperature_2;
	//---------------------------------------

	size_t* left_index = (size_t*) malloc(sizeof(size_t) * size + 1);
	left_index[0] = 1;
	left_index[size] = M - 1;

	size_t i;
	for(i = 1; i < size; i++) 
	{
		left_index[i] = left_index[i - 1] + (M / size) + ((i - 1) < ((M % size) - 2));
	}

	double start = omp_get_wtime();

	// Задаем кол-во процессов
	omp_set_num_threads(size);

	for (n = 0; n < N; n++) 
	{

		#pragma omp parallel
        {
            int id_ = omp_get_thread_num(); // Решил пойти по пути наименьшего сопротивления, чтобы не менять целиком код, поэтому не использовал omp_init/set_lock, а явно разбивать по процессам

			for (m = left_index[id_]; m < left_index[id_ + 1]; m++) 
			{
				u1[m] = u0[m] + 0.3  * (u0[m - 1] - 2.0 * u0[m] + u0[m + 1]);
			}
        }

		double* t = u0;
		u0 = u1;
		u1 = t;
	}
	
    double end = omp_get_wtime();


    
	// Вывод на экран
	printf("%d %.16g\n", size, end - start);//(число процессов, время в секундах)

    
	free(u0);
	free(u1);

	return 0;
}