#include <stdio.h>
#include <omp.h>

int main(int argc, char *argv[])
{
    int n=1;
    printf("serial (1): %d\n", n);
    #pragma omp parallel private(n)
    {
        printf("parallel (1): %d\n", n);
        n = omp_get_thread_num();
        printf("parallel (2): %d\n", n);
    }
    printf("serial (2): %d\n", n);
}
