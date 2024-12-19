#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

double sec(double x)     //функция секанса для производной тангенса
{
	return 1.0 / cos(x);
}

double EPSILON()
{
	double eps = 1.0;
	while (1.0 + eps / 2.0 > 1.0)
	{
		eps = eps / 2.0;
	}
	return eps;
}

double u(double x)
{
	return pow(10, 5) * sin(x); //sin(x); pow(10, 5) * sin(x) ; tan(x);
}

double d3_u(double x)
{
	return -pow(10, 5) * cos(x); //-4 * pow(sec(x), 2) + 6 * pow(sec(x), 4);-pow(10, 5) * cos(x); -cos(x);
}

int main(void)
{
	FILE* file;

	file = fopen("output.txt", "w");

	if (!file)
	{
		printf("Error opening file for writing. \n");
		return 0;
	}

	int    k            = 3;
	double epsilon_mach = EPSILON();

	double x0    = 1.59;
	double u0    = u(x0);
	double d3_u0 = d3_u(x0);

	double d3_u_h;
	double h = pow(10.0, double(k)) * pow(epsilon_mach, 0.25);

	for (int i = k; i >= -k; h /= 10.0, i--)
	{
		d3_u_h = (u(x0 + 2 * h) - 3 * u(x0 + h) + 3 * u0 - u(x0 - h)) / (pow(h, 3));
		
		printf("h = %e, error = %e\n", h, fabs((d3_u_h - d3_u0) / d3_u0));
		fprintf(file,"h = %e, error = %e\n", h, fabs((d3_u_h - d3_u0) / d3_u0));
	}
	
	return 0;
}