#include <stdio.h>
#include <stdlib.h>
#define L 1.0
//T1 -------- T2

#define T1 1.0
#define T2 2.0
//Скорость диффузии = 1


int main(int argc, char * argv[])
{
    if (argc != 3) {
        printf("Usage: %s T M\n", argv[0]);
        return 1;
    }
    double T = atof(argv[1]);
    int M = atoi(argv[2]);
    double h = L / M;
    double tau = 0.3 * h * h;
    int N = T / tau;
    double * u0 = (double *)malloc(sizeof(double)*M);
    double * u1 = (double *)malloc(sizeof(double)*M);
    int i;


    //начальные условия
    for (i = 0; i < M; i++)
    {
        u0[i] = u1[i] + 0.0;
    }


    //граничные условия
    u0[0] = u1[0] = T1;
    u0[M-1] = u1[M-1] = T2;
    int m;
    int n;

    //интегрируем по времени
    for (n = 0; n < N; n++)
    {
        for (m = 1; m < M-1; m++)
        {
            u1[m] = u0[m] + tau / h / h * (u0[m-1] - 2.0*u0[m] + u0[m+1]);
        }
        double *t = u0;
        u0 = u1;
        u1 = t;
    }

    //вывод на экран
    for (m = 0; m < M; m++)
    {
        printf("%lf %lf\n", m*h, u1[i]);
    }
    free(u0);
    free(u1);
    return 0;
}
