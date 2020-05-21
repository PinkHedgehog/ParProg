//декомпозиция по номерам потоков
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>


// Определение функции
double F(double x) {
    return sqrt(4 - x*x);
}

// Формула Котеса рассчета определенного интеграла для равномерной сетки
double Int(size_t left_index, size_t right_index, double h) {
    double I = (F(right_index * h) + F(left_index * h)) / 2;

    for(size_t i = left_index + 1; i < right_index; i++) {
        I += F(i * h);
    }

    return I * h;
}

int main(int argc, char **argv) {
	if (argc != 3) {
		printf("Usage: %s N size\n", argv[0]);
		return 1;
	}
    // Количество шагов
    size_t N = atoi(argv[1]);
    // Запрошенное кол-во процессов
    int size = atoi(argv[2]);
 

    // Задаем границы интегрирования
    double a = 0, b = 2;
    // Задаем мелкость разбиения отрезка
    double h = (b - a) / N;
    double result = 0.0;

        // Начинаем отсчет времени
        double start = omp_get_wtime();

        // Устанавливаем требуемое кол-во процессов
        omp_set_num_threads(size);
        #pragma omp parallel
        {
            // Устанавливаем ранг процесса
            int rank = omp_get_thread_num();
            
            // Передаем каждому процессу "свои" индексы интегрирования
            size_t left_index = rank * (N / size);
            size_t right_index = (rank != size - 1) ? (rank + 1) * (N / size) : N;
            double integral = Int(left_index, right_index, h);

            // Определяем интеграл на заданном интервале
            #pragma omp critical
            {
                result += integral;
            }
        }
  

    // Вывод кол-ва процессов, используемого программой, и усредненное время работы
    printf("%d %lf %lf\n", size, result, omp_get_wtime() - start);

    return 0;
}
