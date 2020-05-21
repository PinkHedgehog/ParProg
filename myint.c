#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

const int G = 3;

double f(double x)
{
    return 4*sqrt(1-x*x);
}

double f2(double x)
{
    return 4/(1+pow(x, 2));
}

int main(int argc, char *argv[])
{
    if (argc != 2)
    {
        printf("Usage: %s <N> \n", argv[0]);
    }
    int p;
    int N = atoi(argv[1]);
    int i;

    int myrank, size;

    MPI_Init(&argc, &argv);
    double start = MPI_Wtime();
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    p = size;
    double l_k = 0;
    double x_i, h;
    h = 1.0 / ((double) N);
    for (i = myrank; i < N; i+= p)
    /*
    формула трапеции для маленького участка: I(a, b, f) ~~ (f(a) + f(b))*(b-a)/2
    */
    {
        x_i = i * h;
        l_k += (f(x_i) + f(x_i + h))*h/2;
    }

    double l = l_k;
    double l_i;
    for (i = 1; i < p; i *= 2)
    {
        if ((myrank % (2*i) == 0) && (myrank + i < p))
        {
             MPI_Recv(&l_i, 1, MPI_DOUBLE, myrank+i, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
             l += l_i;
             //printf("Have %lf at %d\n", l);
        }
        else
        {
            if (myrank % (2*i) == i)
            {
                //printf("Sending %lf from %d to %d\n", l_k, myrank, myrank - i);
                MPI_Send(&l, 1, MPI_DOUBLE, myrank - i, 42, MPI_COMM_WORLD);
            }
        }
    }
    if (myrank == 0)
    {
        l -= (f(0) + f(1))*h/2;
        double end = MPI_Wtime();
        printf("Integral: %.12f, Processes: %d, Partitions: %d, Time: %.5f\n", l, p, N, (end-start));
    }

    // double m;
    // if (myrank == 0)
    // {
    //     m = l_k;
    // }
    // else
    // {
    //     m = 0;
    // }
    // if (myrank % G > 0)
    // {
    //     int tag = myrank < G ? 4 : 3;
    //     printf("Sending %f from %d to %d at A\n", l_k, myrank, myrank - myrank % G );
    //     MPI_Send(&l_k, 1, MPI_DOUBLE, myrank - myrank % G, tag, MPI_COMM_WORLD);
    // }
    // else
    // {
    //     if ((myrank % G == 0) && (myrank > 0))
    //     {
    //         double l = 0;
    //         double l_i;
    //         for (i = myrank + 1; i < myrank + G; i++)
    //         {
    //             if (i >= p)
    //             {
    //                 break;
    //             }
    //             printf("myrank = %d, receiving from %d ", myrank, i);
    //             MPI_Recv(&l_i, 1, MPI_DOUBLE, myrank, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    //             printf(" %f\n", l_i);
    //             l += l_i;
    //         }
    //         l += l_k;
    //         //if (myrank > 0)
    //         printf("Sending %f from %d to %d at B\n", l_k, myrank, 0);
    //         MPI_Send(&l, 0, MPI_DOUBLE, 0, 4, MPI_COMM_WORLD);
    //         //l -= (f(0) + f(1))/2;
    //         //l /= ((double) N);
    //         //double end = MPI_Wtime();
    //         //printf("Integral: %.12f, Processes: %d, Partitions: %d, Time: %.5f, Receiving: %.5f\n", l, p, N, (end-start), (end - rcv));
    //     }
    //
    //     // else
    //     // {
    //     //     MPI_Send(&l_k, 1, MPI_DOUBLE, myrank / G, 3, MPI_COMM_WORLD);
    //     // }
    //
    //
    // }
    // if (myrank == 0)
    // {
    //     double l = l_k;
    //     double l_i;
    //     // for (i = 1; i < p && i < G; i++)
    //     // {
    //     //     MPI_Recv(&l_i, 1, MPI_DOUBLE, i, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    //     //     l += l_i;
    //     // }
    //     for (i = 1; i < p; i += 1)
    //     {
    //         // MPI_Recv(&l_i, 1, MPI_DOUBLE, i, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    //         // l += l_i;
    //         MPI_Recv(&l_i, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    //         l += l_i;
    //     }
    //     l -= (f(0) + f(1))/2;
    //     //t += l;
    //     //l /= ((double) N);
    //     double end = MPI_Wtime();
    //     printf("Integral: %.12f, Processes: %d, Partitions: %d, Time: %.5f\n", l, p, N, (end-start));
    // }

    MPI_Finalize();

    return 0;
}
