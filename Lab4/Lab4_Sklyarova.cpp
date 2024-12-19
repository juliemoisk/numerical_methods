#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define N1 2000
#define N2 4000
#define left_end 0.0
#define right_end 1.0

#define lambda_a 1
#define mu_a     0
#define psi_a    0

#define lambda_b 1
#define mu_b     2
#define psi_b    0

double p(double x)
{
	return x;
}

double q(double x)
{
	return -sqrt(x);
}

double g(double x)
{
	return -3 * exp(-x);
}

//-----------------------------------------------------------------------------------------------------------------

double max(int N, double* array)
{
	double max = array[1];

	for (int i = 2; i <= N; i++)
	{
		if (array[i] > max)
		{
			max = array[i];
		}
	}

	return max;
}

double norm(FILE *file, int N, double* u, double* v)
{
	double *difference = (double*)malloc((N + 1) * sizeof(double));

	if (difference == NULL)
	{
		printf("MEMORY ERROR !\n");
		fprintf(file, "MEMORY ERROR !\n");
		return -1.0;
	}

	double norm;

	for (int i = 1; i <= N; i++)
	{
		difference[i] = fabs(u[i] - v[i]);
	}

	norm = max(N, difference);

	free(difference);

	return norm;
}

double norm_1(FILE *file, int N, double* u, double* v)
{
	double* difference = (double*)malloc((N + 1) * sizeof(double));

	if (difference == NULL)
	{
		printf("MEMORY ERROR !\n");
		fprintf(file, "MEMORY ERROR !\n");
		return -1.0;
	}

	double norm_1 = 0.0;

	for (int i = 1; i <= N; i++)
	{
		difference[i] = fabs(u[i] - v[i]);
		norm_1 += difference[i];
	}

	free(difference);

	return norm_1;
}

double norm_2(FILE *file, int N, double* u, double* v)
{
	double* difference = (double*)malloc((N + 1) * sizeof(double));

	if (difference == NULL)
	{
		printf("MEMORY ERROR !\n");
		fprintf(file, "MEMORY ERROR !\n");
		return -1;
	}

	double sum = 0.0;
	double norm_2;

	for (int i = 1; i <= N; i++)
	{
		difference[i] = pow(fabs(u[i] - v[i]), 2);
		sum += difference[i];
	}

	norm_2 = sqrt(sum);

	free(difference);

	return norm_2;
}

//-----------------------------------------------------------------------------------------------------------------

void sweep_method(FILE *file, int N, double* a, double* b, double* c, double* f, double* u)
{
	double* alpha = (double*)malloc((N + 1) * sizeof(double));
	double* beta =  (double*)malloc((N + 2) * sizeof(double));

	if (alpha == NULL || beta == NULL)
	{
		printf("MEMORY ERROR !\n");
		fprintf(file, "MEMORY ERROR !\n");
		return;
	}

	alpha[1] = b[0] / c[0];

	for (int i = 1; i <= N - 1; i++)
	{
		alpha[i + 1] = b[i] / (c[i] - a[i] * alpha[i]);
	}

	beta[1] = f[0] / c[0];

	for (int i = 1; i <= N; i++)
	{
		beta[i + 1] = (f[i] + a[i] * beta[i]) / (c[i] - a[i] * alpha[i]);
	}

	u[N] = beta[N + 1];

	for (int i = N; i >= 1; i--)
	{
		u[i - 1] = alpha[i] * u[i] + beta[i];
	}

	free(alpha);
	free(beta);
}

void difference_method(int N, double h, double* a, double* b, double* c, double* f, double* p_i, double* q_i, double* g_i)
{
	c[1] = 2 / pow(h, 2) - q_i[1];
	b[1] = 1 / pow(h, 2) + p_i[1] / (2 * h);

	f[1] = (1 / pow(h, 2) - p_i[1] / (2 * h)) * psi_a / lambda_a - g_i[1];

	for (int i = 2; i <= N - 1; i++)
	{
		a[i] = 1 / pow(h, 2) - p_i[i] / (2 * h);
		c[i] = 2 / pow(h, 2) - q_i[i];
		b[i] = 1 / pow(h, 2) + p_i[i] / (2 * h);
		f[i] = -g_i[i];
	}

	a[N] = mu_b / pow(h, 2) - (mu_b * p_i[N - 1]) / h + (mu_b * q_i[N - 1]) / 2;
	c[N] = mu_b / pow(h, 2) - (mu_b * p_i[N - 1] - lambda_b) / h - (lambda_b * p_i[N - 1]) / 2;
	f[N] = (1 / h - p_i[N - 1] / 2) * psi_b - mu_b / 2 * g_i[N - 1];

}

void knot(int N, double h, double* a, double* b, double* c, double* f, double* p_i, double* q_i, double* g_i)
{
	c[1] = 2 / pow(h, 2) - q_i[1];
	b[1] = 1 / pow(h, 2) + p_i[1] / (2 * h);

	f[1] = (1 / pow(h, 2) - p_i[1] / (2 * h)) * psi_a / lambda_a - g_i[1];

	for (int i = 2; i <= N - 1; i++)
	{
		a[i] = 1 / pow(h, 2) - p_i[i] / (2 * h);
		c[i] = 2 / pow(h, 2) - q_i[i];
		b[i] = 1 / pow(h, 2) + p_i[i] / (2 * h);
		f[i] = -g_i[i];
	}

	a[N] = mu_b / pow(h, 2);
	c[N] = mu_b / pow(h, 2) + lambda_b / h + (lambda_b * p_i[N] - mu_b * q_i[N]) / 2;
	f[N] = (1 / h + p_i[N] / 2) * psi_b - mu_b / 2 * g_i[N];
}

void u_and_v(FILE *file, FILE *file_u, int N, double* u, double* v)
{
	double h = (right_end - left_end) / N;

	printf("N = %d, h = %f\n\n", N, h);
	fprintf(file, "N = %d, h = %f\n\n", N, h);

	double* a = (double*)malloc((N + 1) * sizeof(double));
	double* b = (double*)malloc((N + 1) * sizeof(double));
	double* c = (double*)malloc((N + 1) * sizeof(double));
	double* f = (double*)malloc((N + 1) * sizeof(double));

	double* x   = (double*)malloc((N + 1) * sizeof(double));
	double* p_i = (double*)malloc((N + 1) * sizeof(double));
	double* q_i = (double*)malloc((N + 1) * sizeof(double));
	double* g_i = (double*)malloc((N + 1) * sizeof(double));

	if (a == NULL || b == NULL || c == NULL || f == NULL || x == NULL
		|| p_i == NULL || q_i == NULL || g_i == NULL)
	{
		printf("MEMORY ERROR ! \n");
		fprintf(file, "MEMORY ERROR ! \n");
		free(a); free(b); free(c); free(f); free(x); free(p_i); free(q_i); free(g_i);
		return;
	}

	for (int i = 0; i <= N; i++)
	{
		x[i] = left_end + h * i;
		p_i[i] = p(x[i]);
		q_i[i] = q(x[i]);
		g_i[i] = g(x[i]);
	}

	difference_method (N, h, a, b, c, f, p_i, q_i, g_i);
	sweep_method      (file, N, a, b, c, f, u);

	for (int i = 1; i <= N; i++)
	{
		fprintf(file_u, "%f, %f\n", x[i], u[i]);
	}

	knot         (N, h, a, b, c, f, p_i, q_i, g_i);
	sweep_method (file, N, a, b, c, f, v);

	free(a); free(b); free(c); free(f); free(x); free(p_i); free(q_i); free(g_i);
}
//-----------------------------------------------------------------------------------------------------------------

int main(void)
{
	FILE* file;

	fopen_s(&file, "output.txt", "w");

	if (!file)
	{
		printf("Error opening file for writing. \n");
		return 0;
	}
	FILE* file_u1;

	fopen_s(&file_u1, "u1_output.txt", "w");

	if (!file_u1)
	{
		printf("Error opening file for writing. \n");
		return 0;
	}

	FILE* file_u2;

	fopen_s(&file_u2, "u2_output.txt", "w");

	if (!file_u2)
	{
		printf("Error opening file for writing. \n");
		return 0;
	}

	double u1[N1 + 1], v1[N1 + 1], u2[N2 + 1], v2[N2 + 1];
	
	u1[0] = 0.0;
	u2[0] = 0.0;
	v1[0] = 0.0;
	v2[0] = 0.0;

	double norm_n1, norm_1_n1, norm_2_n1;
	double norm_n2, norm_1_n2, norm_2_n2;

	double check, check_1, check_2;

	u_and_v(file, file_u1, N1, u1, v1);
	u_and_v(file, file_u2, N2, u2, v2);

	norm_n1   = norm  (file, N1, u1, v1);
	norm_1_n1 = norm_1(file, N1, u1, v1) / N1;
	norm_2_n1 = norm_2(file, N1, u1, v1) / N1;

	norm_n2   = norm  (file, N2, u2, v2);
	norm_1_n2 = norm_1(file, N2, u2, v2) / N2;
	norm_2_n2 = norm_2(file, N2, u2, v2) / N2;


	printf("N1 = %d: ||u1 - v1|| = %.11f, 1 / N1 * ||u1 - v1||_1 = %.11f, 1 / N1 *||u1 - v1||_2 = %.11f\n", N1, 
		norm_n1, norm_1_n1, norm_2_n1);

	fprintf(file, "N1 = %d: ||u1 - v1|| = %.11f, 1 / N1 * ||u1 - v1||_1 = %.11f, 1 / N1 *||u1 - v1||_2 = %.11f\n", N1,
		norm_n1, norm_1_n1, norm_2_n1);

	printf("N2 = %d: ||u2 - v2|| = %.11f, 1 / N2 * ||u2 - v2||_1 = %.11f, 1 / N2 *||u2 - v2||_2 = %.11f\n", N2,
		norm_n2, norm_1_n2, norm_2_n2);

	fprintf(file, "N2 = %d: ||u2 - v2|| = %.11f, 1 / N2 * ||u2 - v2||_1 = %.11f, 1 / N2 *||u2 - v2||_2 = %.11f\n", N2,
		norm_n2, norm_1_n2, norm_2_n2);

	check   = norm_n1   / norm_n2;
	check_1 = norm_1_n1 / norm_1_n2;
	check_2 = norm_2_n1 / norm_2_n2;

	printf("\ncheck = %f, ckeck_1 = %f, check_2 = %f", check, check_1, check_2);
	fprintf(file,"\ncheck = %f, ckeck_1 = %f, check_2 = %f", check, check_1, check_2);

	return 0;
}
