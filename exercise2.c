#include <stdio.h>
#include <math.h>
#include <stdlib.h>

//定义g(t,x),由题目给出
double g_function(double t, double x) 
{
    return t + sin(x * t * t);
}

//函数e^(x+y)
double exp_xy(double x, double y) 
{
    return exp(x + y);
}

//复化梯形积分
double trapezoidal(double a, double b, int n, double x) 
{
    double h = (b - a) / n;
    double sum = 0.5 * (g_function(a, x) + g_function(b, x));
    for (int k = 1; k < n; k++) 
    {
        double t = a + k * h;
        sum += g_function(t, x);
    }
    return h * sum;
}

//复化Simpson积分
double simpson(double a, double b, int n, double x) 
{
    double h = (b - a) / n;
    double sum = g_function(a, x) + g_function(b, x);
    for (int k = 1; k < n; k += 2) 
    {
        double t = a + k * h;
        sum += 4 * g_function(t, x);
    }
    for (int k = 2; k < n; k += 2) 
    {
        double t = a + k * h;
        sum += 2 * g_function(t, x);
    }
    return h * sum / 3;
}

//计算收敛阶
double calculate_order(double a, double b, int n, double x, int method) 
{
    double I_n, I_2n, I_4n;
    if (method == 0) 
    {
        I_n = trapezoidal(a, b, n, x);
        I_2n = trapezoidal(a, b, 2 * n, x);
        I_4n = trapezoidal(a, b, 4 * n, x);
    } 
    else 
    {
        if (n % 2 != 0) n++;
        I_n = simpson(a, b, n, x);
        I_2n = simpson(a, b, 2 * n, x);
        I_4n = simpson(a, b, 4 * n, x);
    }
    double cal1 = fabs(I_n - I_2n);
    double cal2 = fabs(I_2n - I_4n);
    return log2(cal1 / cal2);  // 按照相对误差公式计算
}

//Romberg积分，同时控制精度
double romberg(double a, double b, double x, double eps) 
{
    int max = 20; 
    double **R = (double **)malloc(max * sizeof(double *));
    for (int i = 0; i < max; i++) 
    {
        R[i] = (double *)malloc(max * sizeof(double));
    }

    R[0][0] = 0.5 * (b - a) * (g_function(a, x) + g_function(b, x));

    for (int k = 1; k < max; k++) 
    {
        int n = 1 << k;
        double h = (b - a) / n;
        double sum = 0;
        for (int i = 1; i < n; i += 2) 
        { 
            double t = a + i * h;
            sum += g_function(t, x);
        }
        R[k][0] = 0.5 * R[k-1][0] + h * sum;

        //Richardson外推
        for (int m = 1; m <= k; m++) 
        {
            double factor = 1.0 / (pow(4, m) - 1);
            R[k][m] = R[k][m-1] + factor * (R[k][m-1] - R[k-1][m-1]);
        }

        if (fabs(R[k][k] - R[k-1][k-1]) < eps) 
        {
            double result = R[k][k];
            for (int i = 0; i < max; i++) 
            {
                free(R[i]);
            }
            free(R);
            return result;
        }
    }

    double result = R[max-1][max-1];
    for (int i = 0; i < max; i++) 
    {
        free(R[i]);
    }
    free(R);
    return result;
}

//二重复化Simpson积分
double double_simpson(double x_a, double x_b, int n_x, double y_a, double y_b, int n_y, double (*f)(double, double)) 
{
    if (n_x % 2 != 0) n_x++;
    if (n_y % 2 != 0) n_y++;
    double h_x = (x_b - x_a) / (double)n_x;  
    double h_y = (y_b - y_a) / (double)n_y;  
    double total = 0.0;

    for (int i = 0; i <= n_x; i++) 
    {
        double x = x_a + (double)i * h_x;  
        double w_x = 0.0;
        if (i == 0 || i == n_x) 
        {
            w_x = 1.0;
        } 
        else if (i % 2 == 1) 
        {
            w_x = 4.0;
        } 
        else 
        {
            w_x = 2.0;
        }

        double sum_y = 0.0;
        for (int j = 0; j <= n_y; j++)
        {
            double y = y_a + (double)j * h_y;  
            double w_y = 0.0;
            if (j == 0 || j == n_y) 
            {
                w_y = 1.0;
            } 
            else if (j % 2 == 1) 
            {
                w_y = 4.0;
            } 
            else 
            {
                w_y = 2.0;
            }
            sum_y += w_y * f(x, y);
        }

        total += w_x * sum_y;
    }
    return total * h_x * h_y / 9.0;
}

//计算内层积分
double f_fun(double x) 
{
    return simpson(0, 1, 1000, x); 
}

//计算外层积分
double simpson_f(double a, double b, int n) 
{
    if (n % 2 != 0) n++;
    double h = (b - a) / n;
    double sum = f_fun(a) + f_fun(b);
    
    for (int k = 1; k < n; k += 2) 
    {
        double x = a + k * h;
        sum += 4 * f_fun(x);
    }
    for (int k = 2; k < n; k += 2) 
    {
        double x = a + k * h;
        sum += 2 * f_fun(x);
    }
    return h * sum / 3;
}


int main() 
{
    //第一题
    printf("=== Task 1 ===\n");
    double a = 0, b = 1, x1 = 1;
    int n_list[] = {10, 20, 40, 80, 160, 320};
    int n_count = sizeof(n_list) / sizeof(n_list[0]);
    for (int i = 0; i < n_count; i++) 
    {
        int n = n_list[i];
        double integral = trapezoidal(a, b, n, x1);
        double order = calculate_order(a, b, n, x1, 0);
        printf("n = %d, Integral = %.12f, Convergence Order = %.4f\n", n, integral, order);
    }

    //第二题
    printf("\n=== Task 2 ===\n");
    for (int i = 0; i < n_count; i++) 
    {
        int n = n_list[i];
        double integral = simpson(a, b, n, x1);
        double order = calculate_order(a, b, n, x1, 1);
        printf("n = %d, Integral = %.12f, Convergence Order = %.4f\n", n, integral, order);
    }

    //第四题
    printf("\n=== Task 4 ===\n");
    double x_list[] = {100, 1000, 10000};
    int x_count = sizeof(x_list) / sizeof(x_list[0]);
    int n_explore[] = {100, 500, 1000, 2000, 5000, 10000};
    int n_explore_count = sizeof(n_explore) / sizeof(n_explore[0]);
    for (int i = 0; i < x_count; i++) 
    {
        double x = x_list[i];
        printf("x = %.0f:\n", x);
        for (int j = 0; j < n_explore_count; j++) 
        {
            int n = n_explore[j];
            double order_trap = calculate_order(a, b, n, x, 0);
            double order_simp = calculate_order(a, b, n, x, 1);
            printf("  n = %d, Trapezoid Order = %.4f, Simpson Order = %.4f\n", n, order_trap, order_simp);
        }
    }

    //第五题
    printf("\n=== Task 5 ===\n");
    double eps = 1e-12;
    double x5_1 = 1, x5_2 = 100;
    double romberg_1 = romberg(a, b, x5_1, eps);
    double romberg_2 = romberg(a, b, x5_2, eps);
    printf("x=1, Romberg Integral = %.12f\n", romberg_1);
    printf("x=100, Romberg Integral = %.12f\n", romberg_2);

    //第六题1
    printf("\n=== Task 6.1 ===\n");
    double x_a1 = 0.0, x_b1 = 3.0, y_a1 = 0.0, y_b1 = 1.0;
    int n_x1_list[] = {6, 12, 24, 48}; 
    int n_y1_list[] = {2, 4, 8, 16};  
    int n_x1_count = sizeof(n_x1_list)/sizeof(n_x1_list[0]);
    double theoretical = 32.794331281498;
    printf("6.1: Theoretical = %.12f\n", theoretical);
    for (int i = 0; i < n_x1_count; i++) 
    {
        int n_x = n_x1_list[i];
        for (int j = 0; j < n_x1_count; j++) 
        {
            int n_y = n_y1_list[j];
            double integral = double_simpson(x_a1, x_b1, n_x, y_a1, y_b1, n_y, exp_xy);
            double error = fabs(theoretical - (integral));
            printf("n_x=%d, n_y=%d, Calculated Integral = %.12f, Error = %.12f\n", n_x, n_y, integral, error);
        }
    }

    //第六题2
    printf("\n=== Task 6.2 ===\n");
    double x_a2 = 0, x_b2 = 10;
    int n_outer_list[] = {2, 4, 8, 16,100,1000};
    int n_outer_count = sizeof(n_outer_list)/sizeof(n_outer_list[0]);
    for (int i = 0; i < n_outer_count; i++) 
    {
        int n_outer = n_outer_list[i];
        double integral = simpson_f(x_a2, x_b2, n_outer);
        printf("Outer n=%d, Integral = %.12f\n", n_outer, -integral);
    }
    return 0;
}