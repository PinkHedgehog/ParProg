#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <omp.h>


double F(double x) 
{
    return sqrt(4 - x*x);
}


double Count_Int(size_t left_index, size_t right_index, double h) 
{
    double I = (F(right_index * h) + F(left_index * h)) / 2;
    for(size_t i = left_index + 1; i < right_index; i++) 
    {
        I += F(i * h);
    }
    return I * h;
}

int main(int argc, char **argv) {

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
    double result = 0;


    omp_set_num_threads(size);

    double start = omp_get_wtime();

    #pragma omp parallel
    {
        int rank = omp_get_thread_num();

        size_t left_index = rank * (N / size);
        size_t right_index = (rank != size - 1) ? (rank + 1) * (N / size) : N;
        double integral = Count_Int(left_index, right_index, h);

            result += integral;
    }

    double end = omp_get_wtime();

    printf("%d %lf %.16g\n", size, result, end - start);

	return EXIT_SUCCESS;
}
