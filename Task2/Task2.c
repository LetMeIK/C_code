#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void real(int k, int n, double A, double *y);
void real(int k, int n, double A, double *y)
{
    double h = pow(10, -n);

        if(k == 1)
        {
            y[0] = 1;

            for(int i = 1; i <= pow(10, n); i++)
                y[i] = y[i-1] - A*h*y[i-1];
        }

        if(k == 2)
        {
            y[0] = 1;

            for(int i = 1; i <= pow(10, n); i++)
                y[i] = y[i-1] / (1.0 + A*h);
        }

        if(k == 3)
        {
            y[0] = 1;

            for(int i = 1; i <= pow(10, n); i++)
                y[i] = y[i-1] * (2.0-A*h) / (2.0 + A*h);
        }

        if(k == 4)
        {
            y[0] = 1;
            y[1] = 1 - A*h;

            for(int i = 2; i <= pow(10, n); i++)
                y[i] = y[i-2] - 2.0*A*h*y[i-1];
        }

        if(k == 5)
        {
            y[0] = 1;
            y[1] = 1 - A*h;

            for(int i = 2; i <= pow(10, n); i++)
                y[i] = (2.0*y[i-1] - 0.5*y[i-2]) / (1.5 + A*h);
        }

        if(k == 6)
        {
            y[0] = 1;
            y[1] = 1 - A*h;

            for(int i = 2; i <= pow(10, n); i++)
                y[i] = 2.0*((A*h-1.5)*y[i-2] + 2.0*y[i-1]);
        }
}


int main(void)
{
    int k, n, p;
    double A, *y;
    printf("Input scheme id: ");
    scanf("%d", &k);
    printf("\nInput n: ");
    scanf("%d", &n);
    printf("\nInput A: ");
    scanf("%lf", &A);

    p = pow(10, n);
    
    y = (double*)malloc((p+1)*sizeof(double));
    real(k, n, A, y);

    p = pow(10, n);

    printf("\n||y(x_n) - y_n|| = %.4e\n", fabs(exp(-A) - y[p]));
    
    return 0;
}
