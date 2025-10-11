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

void generateEquidistantNodes(double* xs, int n) 
{
    for (int i = 0; i < n; i++) 
    {
        xs[i] = -1 + 2.0 * i / n;
    }
}

// 分片线性插值
double piecewiseLinearInterpolate(double x, double* xs, double* ys, int n) 
{
    int i = 0;
    while (i < n - 1 && x > xs[i + 1]) 
    {
        i++;
    }

    double x0 = xs[i];
    double x1 = xs[i + 1];
    double y0 = ys[i];
    double y1 = ys[i + 1];

    return y0 + (y1 - y0) * (x - x0) / (x1 - x0);
}

// 分片 Hermite 插值
double piecewiseHermiteInterpolate(double x, double* xs, double* ys, double* ys_deriv, int n) 
{
    int i = 0;
    while (i < n - 1 && x > xs[i + 1]) 
    {
        i++;
    }

    double x0 = xs[i];
    double x1 = xs[i + 1];
    double y0 = ys[i];
    double y1 = ys[i + 1];
    double d0 = ys_deriv[i];
    double d1 = ys_deriv[i + 1];
    double t = (x - x0) / (x1 - x0);
    double t2 = t * t;
    double t3 = t2 * t;
    double h00 = 2 * t3 - 3 * t2 + 1;
    double h10 = t3 - 2 * t2 + t;
    double h01 = -2 * t3 + 3 * t2;
    double h11 = t3 - t2;
    return h00 * y0 + h10 * (x1 - x0) * d0 + h01 * y1 + h11 * (x1 - x0) * d1;
}

int main() 
{
    int n_values[] = {5, 10, 20, 40};
    int num_n = sizeof(n_values) / sizeof(n_values[0]);
    int N = 200;  // 用于计算误差的稠密节点数

    printf("=== 分片线性插值误差（等距节点）===\n");
    printf("n\tMax Error\n");

    for (int k = 0; k < num_n; k++) 
    {
        int n = n_values[k];
        double* xs = (double*)malloc(n * sizeof(double));
        double* ys = (double*)malloc(n * sizeof(double));

        // 生成等距节点
        generateEquidistantNodes(xs, n);

        // 计算节点处的函数值
        for (int i = 0; i < n; i++) 
        {
            ys[i] = f(xs[i]);
        }

        double max_error = 0.0;
        for (int j = 0; j <= N; j++) 
        {
            double xj = -1 + 2.0 * j / N;
            double y_interp = piecewiseLinearInterpolate(xj, xs, ys, n);
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
    }

    printf("\n=== 分片 Hermite 插值误差（等距节点）===\n");
    printf("n\tMax Error\n");

    for (int k = 0; k < num_n; k++) 
    {
        int n = n_values[k];
        double* xs = (double*)malloc(n * sizeof(double));
        double* ys = (double*)malloc(n * sizeof(double));
        double* ys_deriv = (double*)malloc(n * sizeof(double));

        // 生成等距节点
        generateEquidistantNodes(xs, n);

        // 计算节点处的函数值和导数值
        for (int i = 0; i < n; i++) 
        {
            ys[i] = f(xs[i]);
            ys_deriv[i] = f_derivative(xs[i]);
        }

        double max_error = 0.0;
        for (int j = 0; j <= N; j++) 
        {
            double xj = -1 + 2.0 * j / N;
            double y_interp = piecewiseHermiteInterpolate(xj, xs, ys, ys_deriv, n);
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