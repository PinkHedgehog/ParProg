#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <omp.h>


double F(double x) 
{
    return sqrt(4 - x*x);
}


int main(int argc, char **argv) 
{

    size_t N; // число шагов
    int size; // Число процессов

    if (argc > 1) 
    {
        N = atoll(argv[1]);
		if (argc > 2) 
        {
            size = atoi(argv[2]);
        }
    }


    double a = 0, b = 2;
    double h = (b - a) / N;
    double result = (F(a) + F(b)) / 2;//Будем считать методом трапеций

    omp_set_num_threads(size);
    size_t i;
    //--------------------
    double start = omp_get_wtime();
    //--------------------
    #pragma omp parallel for schedule(static, 10000) reduction(+: result)
        for(i = 1; i < N; i++) 
        {
            result += F(i * h);
        }
    //--------------------
    double end = omp_get_wtime();
    //--------------------
    result *= h;

    printf("%d %lf %.16g\n", size, result, end - start);//(число процессов, значение интеграла, время в секундах)
	return 0;
}