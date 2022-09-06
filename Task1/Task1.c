#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define input "input1.txt" //Файл с собств. знач и векторами
#define matinput "matrix1.txt" // Файл с изначальной матрицей

#define N 3 //Задаем размер матрицы

void InputMatrix(int n, double *a, FILE* file, int inputmode); //Ввод матрицы из файла. (inputmode = 0 - n*n, 1 - n*(n+1))
void InputMatrix(int n, double *a, FILE* file, int inputmode)
{
    for(int i = 0; i < n * (n+inputmode); i++)
    {
        fscanf(file, "%lf", &a[i]);
    }
}

void PrintMatrix(int n, double *a, int inputmode);
void PrintMatrix(int n, double *a, int inputmode)
{
    int i, j;
    for(i = 0; i < n; i++)
    {
        for(j = 0; j < n + inputmode; j++)
        {
            if(a[i*(n+inputmode) + j] < 0)
                printf("%.f ", a[i*(n+inputmode) + j]); // Для красоты положительные числа (т.е. без знака -) сдвинем вправо
            else
                printf(" %.f ", a[i*(n+inputmode) + j]);
        }
        printf("\n");
    }
}

static double Norm(int n, double *a); //Считаем норму вектора
static double Norm(int n, double *a)
{
    double norm = 0;
    for(int i = 0; i < n; i++)
    {
        norm += (a[i] * a[i]);
    }
    norm = sqrt(norm);
    return norm;
}

void MaxScalarProd(int n, double *e); //Вычисляем максимальное скалярное произведение всех пар собств. векторов
void MaxScalarProd(int n, double *e)
{
    double prod = 0, maxprod = 0;
    int max_i, max_j;
    int i, j, k;
    for(i = 0; i < n; i++)
    {
        for(j = i + 1; j < n; j++)
        {
            for(k = 0; k < n; k++)
            {
                prod += (e[i*(n+1) + 1 + k] * e[j*(n+1) + 1 + k]);
            }

            if(abs(prod) >= maxprod)
            {
                maxprod = abs(prod);
                max_i = i;
                max_j = j;
            }

        }
    }
    printf("maxprod = %.3f  |  for (e_%d,e_%d)\n\n", maxprod, max_i + 1, max_j + 1);
}

void MaxDiff(int n, double *a, double *e, double *diff); // Вычисляем максимальную норму ||A*e_i/L_i - e_i|| по всем i
void MaxDiff(int n, double *a, double *e, double *diff)
{
    int i, j ,k;
    int max_i;
    double maxdif = 0;

    for(i = 0; i < n; i++)
    {
        for(j = 0; j < n; j++)
        {
            diff[j] = 0;
            for(k = 0; k < n; k++)
            {
                diff[j] += (a[j*n + k] * e[i*(n+1) + 1 + k]);
            }
            diff[j] /= e[i*(n+1)];
            diff[j] -= e[i*(n+1) + 1 + j];
        }
        if(Norm(n, diff) >= maxdif)
        {
            maxdif = Norm(n, diff);
            max_i = i;
        }

    }
    printf("\nmaxdif  = %.3f  |  for i = %d\n",  maxdif, max_i + 1);

}


int main()
{
    int n = N;

    FILE* fin;
    FILE* fmat;
    fin = fopen(input, "r");
    fmat = fopen(matinput, "r");

    double *matrix; // Изначальная матрица
    double *eigen; // Матрица собственных значений и векторов (структура строки - собств.знач и вектор)
    matrix = (double*)malloc(n * n * sizeof(double));
    eigen = (double*)malloc(n * (n+1) * sizeof(double));

    InputMatrix(n, eigen, fin, 1);
    InputMatrix(n, matrix, fmat, 0);
    fclose(fin);
    fclose(fmat);

    //------------------------------------------------------------------------------------------------------

    double *diff;
    diff = (double*)malloc(n * sizeof(double));

    MaxDiff(n, matrix, eigen, diff);
    MaxScalarProd(n, eigen);

    //------------------------------------------------------------------------------------------------------

    printf("\n\nEigen Values and Vectors\n");
    PrintMatrix(n, eigen, 1);

    printf("\n\nInitial Matrix\n");
    PrintMatrix(n, matrix, 0);

    free(diff);
    free(matrix);
    free(eigen);
    return 0;
}

/* Рабочий пример:
     Matrix1          input1
                    Lmb
     1  1  3        -2  1  0 -1 <- L1, e1
     1  5  1         3  1 -1  1 <- L2, e2
     3  1  1         6  1  2  1 <- L3, e3
*/
