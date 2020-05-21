#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <omp.h>


double f(double x) 
{
    return 4*sqrt(1 - x*x);
}


int main(int argc, char **argv) 
{

    size_t N = 50000000; // число шагов
    int size = 1; // Число процессов

    if (argc > 1) 
    {
        N = atoll(argv[1]);
		if (argc > 2) 
        {
            size = atoi(argv[2]);
        }
    }

/*
    формула трапеции для маленького участка: I(a, b, f) ~~ (f(a) + f(b))*(b-a)/2
    
    {
        x_i = i * h;
        l_k += (f(x_i) + f(x_i + h))*h/2;
    }
*/
    //double a = 0, b = 1;
    double h = 1.0 / N;
    double result = 0; //(f(0) + f(1))/2;

    omp_set_num_threads(size);
    size_t i;
    //--------------------
    double start = omp_get_wtime();
    //--------------------
    #pragma omp parallel for schedule(static, 10000) reduction(+: result)
        for(i = 1; i < N; i++) 
        {
            result += (f(i * h) + f(i*h + h))*h/2;
        }
    //--------------------
    result -= (f(0) + f(1))*h/2;;
    double end = omp_get_wtime();
    //--------------------
    //printf("size,result,time\n");
    printf("%d,%.10lf,%.16g\n", size, result, end - start);
	return 0;
}
