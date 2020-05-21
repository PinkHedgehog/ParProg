#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

double k = 1.0;
//double h = 0.02;
double tau  = 0.0002;
double T = 0.1;


double ut(double u, double u1,double u2, double h, double tau)
{
  return u + (k*tau/(h*h))*(u2 - 2*u + u1);
}

double SolutionU(double x, double t, double l)
{
  double sum = 0;
  int i;
  for (i = 0; i < 100; ++i)
  {
    sum += exp(-(M_PI*M_PI*(2*i + 1)*(2*i + 1)*t/(l*l)))/(2*i + 1)*sin(M_PI*(2*i + 1)*x/l);
  }
  return (sum)*4/M_PI;
}

int KurantCondition(double h, double tau)
{
  //printf("tau = %f, %f, answ %d\n", tau, 0.5*h*h/k, tau < 0.5*h*h/k);
  return (tau < 0.5*h*h/k);
}

void InitU(double *u, long currentN, int myrank, int p)
{
  //printf("Init U, curN = %ld\n", currentN);
  long i;
  for (i = 0; i < currentN + 2; ++i)
  {
    //printf("Init U, i = %d\n", i);
    u[i] = 1;
  }
  //printf("Init U, curN = %ld\n", currentN);
  if (myrank == 0)
  {
    u[0] = 0;
  }
  //printf("Init U, curN = %ld\n", currentN);
  if (myrank == p - 1)
  {
    u[currentN + 1] = 0;
  }
  //printf("Init U, curN = %ld\n", currentN);
}

double *unew(double *u, long currentN, double h, double tau)
{
  //printf("Unew!!!!\n");
  long i = 0;
  double *ucopy = (double *) malloc(sizeof(double)*(currentN + 2));
  for (i = 0; i < currentN + 2; ++i)
  {
    ucopy[i] = u[i];
  }
  for (i = 1; i < currentN + 1; ++i)
  {
    u[i] = ut(ucopy[i], ucopy[i-1] ,ucopy[i+1], h, tau);
  }
  free(ucopy);
  return u;
}

void PrintU(long *printnums, int myrank, long N, int p, double *u)
{
  //printf("PrintU\n");
  int i = 0;
  for (i = 0; i < 11; ++i)
  {
    if (printnums[i] * p /N == myrank)
    {
      printf("x = %f u = %f (myrank  = %d, number = %d, number mod  = %d)\n", 0.1*i, u[printnums[i] % (N/p) + 1], myrank, printnums[i], printnums[i] % (N/p) + 1);
      //printf("Printnum %d, N/p %d, myrank %d, printnums[i] mod (N/p) + 1 = %d\n", printnums[i], N/p, myrank, printnums[i] % (N/p) + 1);
    }
  }
}

void split(long N, int p, long *endnum, long *printnums)
{

  //printf("Split");
  int i;
  long k = 0;
  for (i = 0; i < 11; ++i)
  {
    printnums[i] = (long) (i*0.1*N) - 1;
  }
  printnums[0] = 0;
  for (i = 0; i<p; ++i)
  {
    k += N/p;
    endnum[i] = k - 1;
    if(i == p-1)
    {
      endnum[i] = N - 1;
    }
  }

}


int main(int argc, char *argv[])
{
  if (argc != 3)
    return 1;
  int p = atoi(argv[1]), N = atol(argv[2]);
  double h = 1.0/N;
  /*while (KurantCondition(h, tau) != 1)
  {
    tau /= 10;
  }*/
  tau = h*h/4;
  //printf("h = %f, tau = %f", h, tau);
  long NT = T/tau;
  long endnum[p];
  long printnums[11];
  long currentN, endn, begn;
  split(N, p, endnum, printnums);
  //printf("Main h = %f, tau = %f\n", h, tau);
  int i;
  int myrank, size;
  MPI_Status Status;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  if(myrank == 0)
  {
    printf("-------------- p = %d;   N = %d -----------\n", p, N);
    double start = MPI_Wtime();
    for(i = 1; i < p; ++i)
    {

      if (i == p - 1)
      {
        endn = N - 1;
        begn = endnum[i-1] + 1;
      }
      if (i != p - 1)
      {
        endn = endnum[i];
        begn = endnum[i-1] + 1;
      }
      currentN = endn - begn + 1;
      //printf("Before send begin and end\n");
      if (p > 1)
      {
        MPI_Send(&currentN, 1, MPI_LONG, i, 1 + 10*i, MPI_COMM_WORLD);
        MPI_Send(&begn, 1, MPI_LONG, i, 2 + 10*i, MPI_COMM_WORLD);
        MPI_Send(&endn, 1, MPI_LONG, i, 3 + 10*i, MPI_COMM_WORLD);

      }
      //printf("After send begin and end\n");
    }
    currentN = endnum[0];
    endn = endnum[0];
    begn = 0;
    //preprocessing ends
    //begining calculations

    long t = 0;
    double *u = (double *) malloc(sizeof(long)*(currentN + 2));
    InitU(u, currentN, myrank, p);
    //printf("After init o\n");
    //printf("NT = %ld\n", NT);
    for (t = 0; t < NT; ++t)
    {
      //PrintU(printnums, myrank, N, p, u);
      //sleep(1);
      //printf("-----------------t = %d -------------------------\n",t);
      //printf("u[1] = %f, u[2] = %f, u[3] = %f \n", u[1], u[2], u[3]);
      unew(u, currentN, h, tau);
      //printf("After u[-1] = %f, u[-2] = %f, u[-3] = %f \n", u[currentN+1], u[currentN], u[currentN-1]);
      //PrintU(printnums, myrank, N, p, u);
      //sleep(1);
      //get m-1 and n+1 info
      //printf("Get m-1 and n+1\n");
      if (p > 1)
      {
        MPI_Send(&u[currentN], 1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
        //printf("Myrank %d, during Send, Recv, u[currentN] = %f\n", myrank, u[currentN]);
        MPI_Recv(&u[currentN + 1], 1, MPI_DOUBLE, 1, 2, MPI_COMM_WORLD, &Status);

      }
      //printf("Recieve %f", u[currentN + 1]);
      //printf("Myrank %d, Recieve %f, Send  %f\n", myrank, u[currentN+1], u[currentN]);
      //PrintU(printnums, myrank, N, p, u);

    }
    //PrintU(printnums, myrank, N, p, u);
    double end1 = MPI_Wtime();
    double *u1 = (double *) malloc(sizeof(long)*(N + 2));
    InitU(u1, N, myrank, 1);
    for (t = 0; t < NT; ++t)
    {
      unew(u1, N, h, tau);
    }
    //PrintU(printnums, myrank, N, 1, u1);
    double end2 = MPI_Wtime();
    printf("S = %f, time = %f; time (parallel) = %f\n", (end2-end1)/(end1-start), (end2-end1), (end1-start));
    printf("---------------Answer----------\n");
    for (i = 0; i < 11; ++i)
    {
      printf("x = %f, u = %f\n", i*0.1, SolutionU(0.1*i, T, 1));
    }
    free(u);
    free(u1);
  }
  else
  {
    int endn, begn, currentN;
    double time1 = MPI_Wtime();
    MPI_Recv(&currentN, 1, MPI_LONG, 0, 1 + 10*myrank, MPI_COMM_WORLD, &Status);
    MPI_Recv(&begn, 1, MPI_LONG, 0, 2 + 10*myrank, MPI_COMM_WORLD, &Status);
    MPI_Recv(&endn, 1, MPI_LONG, 0, 3 + 10*myrank, MPI_COMM_WORLD, &Status);
    //printf("Myrank = %d, Recieved info\n", myrank);
    //preprocessing ends
    //begining calculations
    long t = 0;
    double *u = (double *) malloc(sizeof(long)*(currentN + 2));
    InitU(u, currentN, myrank, p);
    //printf("After init 1\n");
    //printf("NT = %Ld\n", NT);
    //PrintU(printnums, myrank, N, p, u);
    for (t = 0; t < NT; ++t)
    {
      //PrintU(printnums, myrank, N, p, u);
      //printf("u[-1] = %f, u[-2] = %f, u[-3] = %f \n", u[currentN], u[currentN-1], u[currentN-2]);
      unew(u, currentN, h, tau);
      //PrintU(printnums, myrank, N, p, u);

      //printf("After u[0] = %f, u[1] = %f, u[2] = %f \n", u[0], u[1], u[2]);
      //get m-1 and n+1 info
      //printf("Myrank = %d, get m - 1, n+1, t = %d from %d\n", myrank, t, NT);
      if (myrank %2 == 0)
      {
        //printf("Warning!!!!!!!!!!!!!!\n");
        if(myrank == p - 1)
        {
          MPI_Send(&u[1], 1, MPI_DOUBLE, myrank - 1, 0, MPI_COMM_WORLD);
          MPI_Recv(&u[0], 1, MPI_DOUBLE, myrank - 1, 2, MPI_COMM_WORLD, &Status);
        }
        else
        {
          MPI_Send(&u[currentN], 1, MPI_DOUBLE, myrank + 1, 0, MPI_COMM_WORLD);
          MPI_Send(&u[1], 1, MPI_DOUBLE, myrank - 1, 0, MPI_COMM_WORLD);
          MPI_Recv(&u[currentN + 1], 1, MPI_DOUBLE, myrank + 1, 2, MPI_COMM_WORLD, &Status);
          MPI_Recv(&u[0], 1, MPI_DOUBLE, myrank - 1, 2, MPI_COMM_WORLD, &Status);
        }

      }
      else
      {
        if(myrank == p - 1)
        {
          //printf("Myrank %d, Before Send, Recv\n", myrank);
          MPI_Recv(&u[0], 1, MPI_DOUBLE, myrank - 1, 0, MPI_COMM_WORLD, &Status);
          //printf("Myrank %d, during Send, Recv, u[0] = %f\n", myrank, u[0]);
          MPI_Send(&u[1], 1, MPI_DOUBLE, myrank - 1, 2, MPI_COMM_WORLD);
          //printf("Myrank %d, during Send, Recv, u[N] = %f\n", myrank, u[currentN]);
          //printf("Myrank %d, Recieve %f, Send  %f\n", myrank, u[0], u[currentN]);
        }
        else
        {
          MPI_Recv(&u[0], 1, MPI_DOUBLE, myrank - 1, 0, MPI_COMM_WORLD, &Status);
          MPI_Recv(&u[currentN + 1], 1, MPI_DOUBLE, myrank + 1, 0, MPI_COMM_WORLD, &Status);
          MPI_Send(&u[currentN], 1, MPI_DOUBLE, myrank + 1, 2, MPI_COMM_WORLD);
          MPI_Send(&u[1], 1, MPI_DOUBLE, myrank - 1, 2, MPI_COMM_WORLD);
        }

      }
      //PrintU(printnums, myrank, N, p, u);
    }
    //PrintU(printnums, myrank, N, p, u);
    //printf("u[1] = %f, u[21] = %f, u[41] = %f\n", u[1], u[21], u[41]);
    double time2 = MPI_Wtime();
    free(u);
    //printf("myrank = %d, N = %d, time = %f\n", myrank, N, (time2-time1));
  }
  MPI_Finalize();

  return 0;
}




//python
/*
import matplotlib
import matplotlib.pyplot as plt
plt.plot([1, 2, 3, 4, 5, 6, 7, 8], [0.65, 0.37, 0.31, 0.11, 0.07, 0.06, 0.06, 0.016], label = '1000')
plt.plot([1, 2, 3, 4, 5, 6, 7, 8], [0.99, 1.99, 2.95, 3.90, 2.66, 2.92, 1.48, 1.8], label = '1000000')
plt.plot([1,2,3,4,5,6,7,8], [0.99, 2, 2.98, 3.97, 3.28, 3.51, 3.55, 3.61], label = '100000000')
plt.grid()
plt.title("S(p)")
plt.savefig("hw1.png")



*/
