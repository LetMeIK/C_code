#include <stdio.h>
#include <stdlib.h>
#include <math.h>


double p(double x);
double p(double x)
{
	return 1.0;
}


double y(double x);
double y(double x)
{
	return sin(M_PI * x) + M_PI * x;
}


double f(double x);
double f(double x)
{
	return M_PI * M_PI * sin(M_PI * x) + sin(M_PI * x) + M_PI * x;
}


double EigenVal(int m, double h, int N);
double EigenVal(int m, double h, int N)
{
	double lambda_n;
	lambda_n = 4.0 * N * N * sin(M_PI * (2.0 * m - 1) * h * 0.25) * sin(M_PI * (2.0 * m - 1) * h * 0.25);

	return lambda_n + p(N);
}


double EigenVec(int m, int k, int N);
double EigenVec(int m, int k, int N)
{
	double c;
	c = sqrt(2.0);

	return c * sin(M_PI * (2 * m - 1) * k / (2.0 * N));
}


double ScalarProd(double *y1, double *y2, int N);
double ScalarProd(double *y1, double *y2, int N)
{
	double res = 0;

	for (int i = 1; i < N; i++)
	{
		res += y1[i] * y2[i];
	}
	res += y1[N] * y2[N] * 0.5;
	res /= N;

	return res;
}

void Progon(double *y, int N);
void Progon(double *y, int N)
{
	double *alpha, *beta;
	double h = 1.0 / N;

	double c_1 = 2.0 / (h * h) + p(h);
	double b_1 = 1.0 / (h * h);
	double f_1 = f(h);

	alpha = (double*)malloc((N + 1) * sizeof(double));
	beta  = (double*)malloc((N + 1) * sizeof(double));

	alpha[2] = b_1 / c_1;
	beta[2]  = f_1 / c_1;

	for (int k = 2; k < N; k++)
	{
		double b_k = 1.0 / (h * h);
		double c_k = 2.0 / (h * h) + p(k * h);
		double a_k = 1.0 / (h * h);
		double f_k = f(k * h);

		alpha[k + 1] = b_k / (c_k - alpha[k] * a_k);
		beta[k + 1]  = (a_k * beta[k] + f_k) / (c_k - alpha[k] * a_k);
	}

	double a_N = 2.0 / (h * h);
	double c_N = 2.0 / (h * h) + p(N * h);

	y[N] = (f(N * h) + a_N * beta[N]) / (c_N - a_N * alpha[N]);

	for (int k = N - 1; k >= 1; k--)
		y[k] = alpha[k + 1] * y[k + 1] + beta[k + 1];

	free(alpha);
	free(beta);
}

void Fourier(double *y, int N);
void Fourier(double *y, int N)
{
	double h, c_m, *f_k, *y_m;
	int i, m;

	h = 1.0 / N;
	f_k = (double*)malloc((N + 1) * sizeof(double));
	y_m = (double*)malloc((N + 1) * sizeof(double));
	for (i = 1; i <= N; i++)
	{
		y[i] = 0;
		f_k[i] = f(i * h);
	}
	for (m = 1; m <= N; m++)
	{
		for (int i = 1; i <= N; i++)
			y_m[i] = EigenVec(m, i, N);
		c_m = ScalarProd(f_k, y_m, N) / EigenVal(m, h, N);
		for (i = 1; i <= N; i++)
			y[i] += c_m * y_m[i];
	}

	free(f_k);
	free(y_m);
}

void y_sol(double *Y, int N);
void y_sol(double *Y, int N)
{
	double h = 1.0 / N;
	for (int i = 0; i <= N; i++)
		Y[i] = y(i * h);
}

double Err(double *y_k, double *y, int N);
double Err(double *y_k, double *y, int N)
{
	double res = 0.0;

	for (int i = 1; i <= N - 1; i++)
	{
		res += (y_k[i] - y[i]) * (y_k[i] - y[i]);
	}
	res += 0.5 * (y_k[N] - y[N]) * (y_k[N] - y[N]);
	return sqrt(res / N);
}

int main(void)
{
	double *y_k, *y;
	int N;
	int Method;
	FILE *f;

	f = fopen("result.txt", "w");
	printf("Enter\n1 for Fourier Method\n2 for Run-Through Method\n");
	scanf("%d", &Method);
	printf("Enter N:\n");
	scanf("%d", &N);

	y_k = (double *)malloc((N + 1) * sizeof(double));
	y   = (double *)malloc((N + 1) * sizeof(double));

	if (Method == 1)
		Fourier(y_k, N);
	else
		Progon(y_k, N);

	y_sol(y, N);
	for (int i = 1; i <= N; i++)
		fprintf(f, "%e %e\n", y_k[i], y[i]);

	double error = Err(y_k, y, N);
	printf("Error = %e\n", error);

	free(y); free(y_k);
	fclose(f);

	return 0;
}
