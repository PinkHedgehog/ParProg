#include <stdio.h>

int main()
{
    #ifdef _OPENMP
        printf("OpenMP is supported, version %d\n", _OPENMP);
    #endif
    printf("Serial\n");
#pragma omp parallel num_threads(13)
{
    printf("Parallel\n");
}
    printf("Serial\n");
    return 0;
}
