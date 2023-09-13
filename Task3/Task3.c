#include <stdio.h>
#include <stdlib.h>
#include <math.h>


void fil(double *yk, double xk);
void fil(double *yk, double xk)
{
    yk[0] = cos(xk);
    yk[1] = sin(xk);
}


double *func(double *f, double x, double *y);
double *func(double *f, double x, double *y)
{
    f[0] = -sin(x);
    f[1] =  cos(x);

    return f;
}


double norm(double *a, double *b, int N);
double norm(double *a, double *b, int N)
{
    double s = 0.0;

    for(int i = 0; i < N; i++)
        s += (b[i] - a[i]) * (b[i] - a[i]);

    return sqrt(s);
}


void Step(double *yk, double *tmp, double *yk1, double xk, double h, int N, double *k1, double *k2, double *k3, double *k4, double *f);
void Step(double *yk, double *tmp, double *yk1, double xk, double h, int N, double *k1, double *k2, double *k3, double *k4, double *f)
{
    for(int i = 0; i < N; i++)
    {
        k1[i] = h * func(f, xk, yk)[i];
        tmp[i] = yk[i] + h / 2 * k1[i];
        k2[i] = h * func(f, xk + 0.5 * h, tmp)[i];
        tmp[i] = yk[i] + h / 2 * k2[i];
        k3[i] = h * func(f, xk + 0.5 * h, tmp)[i];
        tmp[i] = yk[i] + h * k1[i];
        k4[i] = h * func(f, xk + h, tmp)[i];
        yk1[i] = (1.0 / 6.0) * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) + yk[i];
    }
}

double RungeErr(double *yk, double *tmp, double *yk1, double *yk21, double *yk22, double xk, double h, int N, double *k1, double *k2, double *k3, double *k4, double *f);
double RungeErr(double *yk, double *tmp, double *yk1, double *yk21, double *yk22, double xk, double h, int N, double *k1, double *k2, double *k3, double *k4, double *f)
{
    Step(yk, tmp, yk1, xk, h, N, k1, k2, k3, k4, f);
    Step(yk, tmp, yk21, xk, h / 2.0, N, k1, k2, k3, k4, f);
    Step(yk21, tmp, yk22, xk + h / 2.0, h / 2.0, N, k1, k2, k3, k4, f);

    return norm(yk1, yk22, N) / (pow(2, 4) - 1);
}

int main(void)
{
    int N = 2, i, cnt;
    double n, h, xk = 0.0, err, e1 = 1e-4, e2 = 1e-3, *yk1, *yk21, *yk22, *yk, *tmp, *y_k_real, *k1, *k2, *k3, *k4, *f;

    scanf("%lf", &n);
    printf("\n");

    h = pow(10.0, -n);

y_k_real  = (double*)malloc(N * sizeof(int));
    yk    = (double*)malloc(N * sizeof(int));
    tmp   = (double*)malloc(N * sizeof(int));
    yk1   = (double*)malloc(N * sizeof(int));
    yk21  = (double*)malloc(N * sizeof(int));
    yk22  = (double*)malloc(N * sizeof(int));
    k1    = (double*)malloc(N * sizeof(int));
    k2    = (double*)malloc(N * sizeof(int));
    k3    = (double*)malloc(N * sizeof(int));
    k4    = (double*)malloc(N * sizeof(int));
    f     = (double*)malloc(N * sizeof(int));

    fil(yk, xk);
    fil(y_k_real, xk);

    printf("                    |      ESTIMATED      |         Real        |            \n");
    printf("   xk    |     h    |   yk_0   |   yk_1   |  y_0(xk) |  y_1(xk) |    error   \n");
    printf("%f | %f | %f | %f | %f | %f | %.3e\n", xk, h,  yk[0], yk[1], y_k_real[0], y_k_real[1], norm(yk, y_k_real, N));


    for(cnt = 0; xk < 1.0; cnt++)
    {
        err = RungeErr(yk, tmp, yk1, yk21, yk22, xk, h, N, k1, k2, k3, k4, f);

        while(err > e2)
        {
            h /= 2.0;
            err = RungeErr(yk, tmp, yk1, yk21, yk22, xk, h, N, k1, k2, k3, k4, f);
        }

        if((xk + h) > 1.0)
        {
            h = 1.0 - xk;
            err = RungeErr(yk, tmp, yk1, yk21, yk22, xk, h, N, k1, k2, k3, k4, f);
        }

        fil(y_k_real, xk + h);
        printf("%f | %f | %f | %f | %f | %f | %.3e\n", xk + h, h, yk22[0], yk22[1], y_k_real[0], y_k_real[1], norm(y_k_real, yk22, N));

        xk += h;
        for(i = 0; i < N; i++)
            yk[i] = yk22[i];

        if(err < e1 && cnt != 0)
            h *= 2;
    }

    printf("\nAccuracy:\ne1 = %.1e \ne2 = %.1e\n\n", e1, e2);
    printf("Dimension of the problem = %d\n", N);
    printf("Number of grid nodes \t = %d\n", cnt);

    free(y_k_real); free(yk); free(tmp); free(yk1); free(yk21); free(yk22); free(k1); free(k2); free(k3); free(k4); free(f);

    return 0;
}
