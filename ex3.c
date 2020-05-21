#include <stdio.h>
#include <omp.h>

int main(int argc, char * argv[])
{
#pragma omp parallel
    {
        printf("Msg 1\n");
#pragma omp barrier
        printf("Msg 2\n");
    }
    return 0;

}
