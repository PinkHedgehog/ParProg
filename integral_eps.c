/*
 * Вычисление интеграла с заданной точностью.
 * Хохлов Николай Игоревич, k_h@inbox.ru, 2014.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

const double A = 0.0001;
const double B = 1;

double function(double x)
{
	double s = sin(1.0 / x ) / x ;
    return s*s;
}

double int_trap(double a, double b, double fa, double fb, double eps)
{
	double I = 0.0;
	double c = (a + b) / 2.0;
	double fc = function(c);
	double sab = (fa + fb) * (b - a) / 2.0;
	double sac = (fa + fc) * (c - a) / 2.0;
	double scb = (fc + fb) * (b - c) / 2.0;
	double sabc = sac + scb;
	if (fabs(sab - sabc) >= eps * fabs(sabc)) {
		I = int_trap(a, c, fa, fc, eps) + int_trap(c, b, fc, fb, eps);
	} else {
		I = sabc;
	}
	return I;
}

int main(int argc, char* argv[])
{
	double eps = 1e-6;
	if (argc > 1) {
	eps = atof(argv[1]);
	}
	printf("%.9lf\n", int_trap(A, B, function(A), function(B), eps));
	return 0;
}
