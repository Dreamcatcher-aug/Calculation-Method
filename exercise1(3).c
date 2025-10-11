#include <stdio.h>
#include <math.h>
#include <stdlib.h>

// 实现拉格朗日插值法
double lagrange(double x, double* xs, double* ys, int n) 
{
    double result = 0.0;
    for(int i = 0; i < n; i++) 
    {
        double L = 1.0;
        for(int j = 0; j < n; j++) 
        {
            if(j != i) 
            {
                L *= (x - xs[j]) / (xs[i] - xs[j]);
            }
        }
        result += L * ys[i];
    }
    return result;
}

// 生成切比雪夫节点
void generateChebyshevNodes(double* xs, int n) 
{
    for (int k = 0; k < n; k++) 
    {
        xs[k] = cos((2 * k + 1) * M_PI / (2 * (n)));
    }
}

// 第一个待验证的函数
double f1(double x) 
{
    return sqrt(1 - x * x);
}

// 第二个待验证的函数
double f2(double x) 
{
    return 1.0 / (1 + 25 * x * x);
}

int main() 
{
    int n_values[] = {5, 10, 20, 40};
    int num_n = sizeof(n_values) / sizeof(n_values[0]);
    int N = 200;  // 按照题目要求使用稠密等距节点，N=200

    // 第一个函数
    printf("=== 函数 f(x) = sqrt(1 - x²)（切比雪夫节点）的插值误差 ===\n");
    printf("n\tMax Error\n");

    for (int k = 0; k < num_n; k++) 
    {
        int n = n_values[k];
        double* xs = (double*)malloc(n * sizeof(double));
        double* ys = (double*)malloc(n * sizeof(double));

        // 生成切比雪夫节点
        generateChebyshevNodes(xs, n);

        for (int i = 0; i < n; i++) 
        {
            ys[i] = f1(xs[i]);       // 计算节点对应的函数值
        }

        double max_error = 0.0;
        for (int j = 0; j <= N; j++) 
        {
            double xj = -1 + 2.0 * j / N;
            double y_interp = lagrange(xj, xs, ys, n);
            double y_true = f1(xj);
            double error = fabs(y_interp - y_true);
            if (error > max_error) 
            {
                max_error = error;
            }
        }

        printf("%d\t%.3e\n", n, max_error);

        // 当n=20时，生成绘图数据
        if (n == 20) 
        {
            FILE* fp = fopen("f1_plot_data.txt", "w");
            for (int j = 0; j <= N; j++) 
            {
                double xj = -1 + 2.0 * j / N;
                double y_interp = lagrange(xj, xs, ys, n);
                double y_true = f1(xj);
                fprintf(fp, "%f %f %f\n", xj, y_true, y_interp);
            }
            fclose(fp);
            printf("已生成 f(x) = sqrt(1 - x²)（切比雪夫节点，n=20）的绘图数据到 f1_plot_data.txt\n");
        }

        free(xs);
        free(ys);
    }

    // 第二个函数
    printf("\n=== 函数 f(x) = 1 / (1 + 25x²)（切比雪夫节点）的插值误差 ===\n");
    printf("n\tMax Error\n");

    for (int k = 0; k < num_n; k++) 
    {
        int n = n_values[k];
        double* xs = (double*)malloc(n * sizeof(double));
        double* ys = (double*)malloc(n * sizeof(double));

        // 生成切比雪夫节点
        generateChebyshevNodes(xs, n);

        for (int i = 0; i < n; i++) 
        {
            ys[i] = f2(xs[i]);
        }

        double max_error = 0.0;
        for (int j = 0; j <= N; j++) 
        {
            double xj = -1 + 2.0 * j / N;
            double y_interp = lagrange(xj, xs, ys, n);
            double y_true = f2(xj);
            double error = fabs(y_interp - y_true);
            if (error > max_error) 
            {
                max_error = error;
            }
        }

        printf("%d\t%.3e\n", n, max_error);

        if (n == 20) 
        {
            FILE* fp = fopen("f2_plot_data.txt", "w");
            for (int j = 0; j <= N; j++) 
            {
                double xj = -1 + 2.0 * j / N;
                double y_interp = lagrange(xj, xs, ys, n);
                double y_true = f2(xj);
                fprintf(fp, "%f %f %f\n", xj, y_true, y_interp);
            }
            fclose(fp);
            printf("已生成 f(x) = 1 / (1 + 25x²)（切比雪夫节点，n=20）的绘图数据到 f2_plot_data.txt\n");
        }

        free(xs);
        free(ys);
    }
    return 0;
}