#include <stdio.h>
#include <math.h>
#include <stdlib.h>

double f(double x) 
{
    return 1.0 / (1 + 25 * x * x);
}

double f_derivative(double x) 
{
    return -50 * x / ((1 + 25 * x * x) * (1 + 25 * x * x));
}

//生成等距节点
void generateEquidistantNodes(double* xs, int n) 
{
    for (int i = 0; i < n; i++) 
    {
        xs[i] = -1 + 2.0 * i / n;
    }
}

//计算Hermite基函数
double hermite_basis(int i, double x, double* xs, int n) 
{
    double result = 1.0;
    double sum = 0.0;
    for (int j = 0; j < n; j++) 
    {
        if (j != i) 
        {
            double term = (x - xs[j]) / (xs[i] - xs[j]);
            result *= term * term;
            sum += 1.0 / (xs[i] - xs[j]);
        }
    }
    return result * (1 - 2 * (x - xs[i]) * sum);
}

//计算Hermite插值
double hermite_interpolate(double x, double* xs, double* ys, double* ys_deriv, int n) 
{
    double result = 0.0;
    for (int i = 0; i < n; i++) 
    {
        double basis = hermite_basis(i, x, xs, n);
        result += ys[i] * basis + (x - xs[i]) * ys_deriv[i] * basis;
    }
    return result;
}

int main()
{
    int n_values[] = {5, 10, 20, 40};
    int num_n = sizeof(n_values) / sizeof(n_values[0]);
    int N = 200;  // 用于计算误差的稠密节点数

    printf("=== 函数 f(x) = 1 / (1 + 25x²) 的 Hermite 插值误差（等距节点）===\n");
    printf("n\tMax Error\n");

    for (int k = 0; k < num_n; k++) 
    {
        int n = n_values[k];
        double* xs = (double*)malloc(n * sizeof(double));
        double* ys = (double*)malloc(n * sizeof(double));
        double* ys_deriv = (double*)malloc(n * sizeof(double));

        generateEquidistantNodes(xs, n);

        for (int i = 0; i < n; i++) 
        {
            ys[i] = f(xs[i]);
            ys_deriv[i] = f_derivative(xs[i]);
        }

        double max_error = 0.0;
        for (int j = 0; j <= N; j++) 
        {
            double xj = -1 + 2.0 * j / N;
            double y_interp = hermite_interpolate(xj, xs, ys, ys_deriv, n);
            double y_true = f(xj);
            double error = fabs(y_interp - y_true);
            if (error > max_error) 
            {
                max_error = error;
            }
        }

        printf("%d\t%.3e\n", n, max_error);

        free(xs);
        free(ys);
        free(ys_deriv);
    }

    return 0;
}