#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define EPSILON    1e-4
#define e          0.00005
#define left_end   0.4
#define right_end  0.6
#define bottom_end 0.6
#define upper_end  0.8

double max(double a, double b) 
{
    return (a > b) ? a : b;
}

//--------------------------------------------------------------------------------------------------------------------

double static f(double x, double y)
{
    return 2 * x * x + y * y - 1;
} 

double static g(double x, double y)
{
    return pow(x, 2. / 3.) - y;
} 

double static phi(double y)
{
    return sqrt(1.0 / 2.0 - y * y / 2.0);
}

double static psi(double x)
{
    return pow(x, 2.0 / 3.0);
}

//--------------------------------------------------------------------------------------------------------------------

double static f_x(double x)
{
    return 4 * x;
}

double static f_y(double y)
{
    return 2 * y;
}

double static g_x(double x)
{
    return 2 * pow(x, -1.0 / 3.0) / 3.0;
}

double static g_y(double y)
{
    return -1;
}

//--------------------------------------------------------------------------------------------------------------------

double static phi_x()
{
    return 0;
}

double static phi_y(double y)
{
    return (-y) / (sqrt(2 * (1 - y * y)));
}

double static psi_x(double x)
{
    return 2 * pow(x, -1.0 / 3.0) / 3.0;
}

double static psi_y()
{
    return 0;
}

//--------------------------------------------------------------------------------------------------------------------

double static f_xx(double x)
{
    return 4;
}

double static f_yy(double y)
{
    return 2;
}

double static g_xx(double x)
{
    return -2 * pow(x, -4.0 / 3.0) / 9.0;
}

double static g_yy(double y)
{
    return 0;
}

double static phi_yy(double y)
{
    return -sqrt(1 - y * y) / (sqrt(2) * (pow(y, 4) - 2 * y * y + 1));
}

double static psi_xx(double x)
{
    return -2 * pow(x, -4.0 / 3.0) / 9.0;
}

//---------------------------------------------------------------------------------------------------------------------

double gradient_ascent (double (*func)(double), double (*derivative_f)(double), double a, double b, FILE *file)
{
    double function_max = -1000;
    double x_max = -1000;

    for (double x = a; x <= b; x += 0.1)
    {
        double x1 = x;
        double x2 = (a + b) / 2;
        double temp = EPSILON + 1;

        while (x1 >= a && x1 <= b && temp >= EPSILON)
        {
            x2 = x1;
            x1 += derivative_f(x1) * 0.001;
            temp = fabs(x1 - x2);

            if (function_max < func(x1) && x1 >= a && x1 <= b)
            {
                function_max = func(x1);
                x_max = x1;
            }
        }
    }
    printf("\tmaximum is in %f and equals %f\n", x_max, function_max);
    fprintf(file, "\tmaximum is in %f and equals %f\n", x_max, function_max);
    return function_max;
}

double gradient_descent(double (*func)(double), double (*derivative_f)(double), double a, double b, FILE *file)
{
    double function_min = 1000;
    double x_min = 1000;

    for (double x = a; x <= b; x += 0.1)
    {
        double x1 = x;
        double x2 = (a + b) / 2;
        double temp = EPSILON + 1;

        while (x1 >= a && x1 <= b && temp >= EPSILON)
        {
            x2 = x1;
            x1 -= derivative_f(x1) * 0.001;
            temp = fabs(x1 - x2);

            if (function_min > func(x1) && x1 >= a && x1 <= b)
            {
                function_min = func(x1);
                x_min = x1;
            }
        }
    }

    printf("\tminimum is in %f and equals %f\n", x_min, function_min);
    fprintf(file,"\tminimum is in %f and equals %f\n", x_min, function_min);

    return function_min;
}

//--------------------------------------------------------------------------------------------------------------------
// оценка нормы матрицы Якоби для (phi(x, y); psi(x, y)), т. е. поиск q

double J_phi_norm(FILE *file) 
{
    printf ("PHI_Y:\n");
    fprintf(file, "PHI_Y:\n");
    double phi_y_max = gradient_ascent (phi_y, phi_yy, bottom_end, upper_end, file);
    double phi_y_min = gradient_descent(phi_y, phi_yy, bottom_end, upper_end, file);

    printf ("\nPSI_X:\n");
    fprintf(file, "\nPSI_X:\n");
    double psi_x_max = gradient_ascent (psi_x, psi_xx, left_end, right_end, file);
    double psi_x_min = gradient_descent(psi_x, psi_xx, left_end, right_end, file);

    double M_phi_y, m_phi_y;
    double M_psi_x, m_psi_x;

    double q;

    if (fabs(phi_y_min) > fabs(phi_y_max))
    {
        M_phi_y = fabs(phi_y_min);
        m_phi_y = fabs(phi_y_max);
    }
    else
    {
        M_phi_y = fabs(phi_y_max);
        m_phi_y = fabs(phi_y_min);
    }

    if (fabs(psi_x_min) > fabs(psi_x_max))
    {
        M_psi_x = fabs(psi_x_min);
        m_psi_x = fabs(psi_x_max);
    }
    else
    {
        M_psi_x = fabs(psi_x_max);
        m_psi_x = fabs(psi_x_min);
    }

    if (M_phi_y > M_psi_x) q = M_phi_y;
    else q = M_psi_x;

    printf ("\n||J|| = %f = q\n", q);
    fprintf(file, "\n||J|| = %f = q\n", q);

    return q;
}

// оценка нормы матрицы Якоби для (f(x, y); g(x, y))

double J_f_norm(FILE *file) 
{
    printf ("\nF_X:\n");
    fprintf(file, "\nF_X:\n");
    double max_f_x = gradient_ascent (f_x, f_xx, left_end, right_end, file);
    double min_f_x = gradient_descent(f_x, f_xx, left_end, right_end, file);

    printf ("\nF_Y:\n");
    fprintf(file, "\nF_Y:\n");
    double max_f_y = gradient_ascent (f_y, f_yy, bottom_end, upper_end, file);
    double min_f_y = gradient_descent(f_y, f_yy, bottom_end, upper_end, file);

    printf ("\nG_X:\n");
    fprintf(file, "\nG_X:\n");
    double max_g_x = gradient_ascent (g_x, g_xx, left_end, right_end, file);
    double min_g_x = gradient_descent(g_x, g_xx, left_end, right_end, file);

    printf ("\nG_Y:\n");
    fprintf(file,"\nG_Y:\n");
    double max_g_y = gradient_ascent (g_y, g_yy, bottom_end, upper_end, file);
    double min_g_y = gradient_descent(g_y, g_yy, bottom_end, upper_end, file);

    double M_f_x, m_f_x;
    double M_f_y, m_f_y;

    double M_g_x, m_g_x;
    double M_g_y, m_g_y;

    double J_f_norm;

    if (fabs(min_f_x) > fabs(max_f_x))
    {
        M_f_x = fabs(min_f_x);
        m_f_x = fabs(max_f_x);
    }
    else
    {
        M_f_x = fabs(max_f_x);
        m_f_x = fabs(min_f_x);
    }

    if (fabs(min_f_y) > fabs(max_f_y))
    {
        M_f_y = fabs(min_f_y);
        m_f_y = fabs(max_f_y);
    }
    else
    {
        M_f_y = fabs(max_f_y);
        m_f_y = fabs(min_f_y);
    }

    if (fabs(min_g_x) > fabs(max_g_x))
    {
        M_g_x = fabs(min_g_x);
        m_g_x = fabs(max_g_x);
    }
    else
    {
        M_g_x = fabs(max_g_x);
        m_g_x = fabs(min_g_x);
    }

    if (fabs(min_g_y) > fabs(max_g_y))
    {
        M_g_y = fabs(min_g_y);
        m_g_y = fabs(max_g_y);
    }
    else
    {
        M_g_y = fabs(max_g_y);
        m_g_y = fabs(min_g_y);
    }

    J_f_norm = max(M_f_x + M_f_y, M_g_x + M_g_y);

    return J_f_norm;
}

//--------------------------------------------------------------------------------------------------------------------
// определитель матрицы Якоби для (f(x, y); g(x, y))

double det_J_f(double x, double y) 
{
    return -4 * x - 2 * y * 2 * pow(x, -1.0 / 3.0) / 3.0;
}

// вычисление минимального значения этого определителя

double det_min() 
{
    return det_J_f(right_end, upper_end);
}

//--------------------------------------------------------------------------------------------------------------------
// вычисление элементов обратной матрицы Якоби

double reverse_J_11(double x, double y)
{
    return g_y(y) / det_J_f(x, y);
}

double reverse_J_12(double x, double y)
{
    return -1 * f_y(y) / det_J_f(x, y);
}

double reverse_J_21(double x, double y)
{
    return -1 * g_x(x) / det_J_f(x, y);
}

double reverse_J_22(double x, double y)
{
    return f_x(x) / det_J_f(x, y);
}

//--------------------------------------------------------------------------------------------------------------------
// вычисление максимумов элементов обратной матрицы Якоби

double reverse_J_11_min(double y)
{
    return g_y(y) / det_min();
}

double reverse_J_12_min(double y)
{
    return -1 * f_y(y) / det_min();
}

double reverse_J_21_min(double x)
{
    return -1 * g_x(x) / det_min();
}

double reverse_J_22_min(double x)
{
    return f_x(x) / det_min();
}

//---------------------------------------------------------------------------------------------------------------------
// оценка нормы обратной матрицы Якоби

double reverse_J_norm()
{
    double fabs_J_11 = fabs(reverse_J_11_min(0.8));
    double fabs_J_12 = fabs(reverse_J_12_min(0.8));
    double fabs_J_21 = fabs(reverse_J_21_min(0.6));
    double fabs_J_22 = fabs(reverse_J_22_min(0.6));

    double reverse_J_norm = max(fabs_J_11 + fabs_J_12, fabs_J_21 + fabs_J_22);

    return reverse_J_norm;
}

double mu_find(FILE * file)
{
    return J_f_norm(file) * reverse_J_norm();
}

//--------------------------------------------------------------------------------------------------------------------

double simple_iteration_method(double q, FILE *file)
{
    int iterations = 0;

    double x0 = left_end;
    double y0 = bottom_end;

    double x, y, nullity;

    while (1)
    {
        x = phi(y0);
        y = psi(x0);

        if (max(fabs(x - x0), fabs(y - y0)) < (1 - q) * EPSILON / q)
        {
            break;
        }
        
        x0 = x;
        y0 = y;
        iterations++;

        if (iterations > 10000)
        {
            printf ("\n Too many iterations !");
            fprintf(file,"\n Too many iterations !");
            return -2;
        }
    }

    nullity = max(fabs(f(x, y)), fabs(g(x, y)));

    printf ("\nSimple iteration method gave the root (x, y) = (%f, %f) in %d iterations and ||f(x, y); g(x, y)|| = %f\n", x, y, iterations + 1, nullity);
    fprintf(file,"\nSimple iteration method gave the root (x, y) = (%f, %f) in %d iterations and ||f(x, y); g(x, y)|| = %f\n", x, y, iterations + 1, nullity);
    return nullity;
}

double seidel_method(FILE *file)
{
    int iterations = 0;

    double x0 = left_end;
    double y0 = bottom_end;

    double x, y, nullity;

    while (1)
    {
        x = phi(y0);
        y = psi(x);

        if (max(fabs(x - x0), fabs(y - y0)) < EPSILON)
        {
            break;
        }

        x0 = x;
        y0 = y;
        iterations++;

        if (iterations > 10000)
        {
            printf ("\n Too many iterations !");
            fprintf(file,"\n Too many iterations !");
            return -2;
        }
    }

    nullity = max(fabs(f(x, y)), fabs(g(x, y)));

    printf ("\nSeidel method gave the root (x, y) = (%f, %f) in %d iterations and ||f(x, y); g(x, y)|| = %f\n", x, y, iterations + 1, nullity);
    fprintf(file, "\nSeidel method gave the root (x, y) = (%f, %f) in %d iterations and ||f(x, y); g(x, y)|| = %f\n", x, y, iterations + 1, nullity);
    return nullity;
}

double newton_method(double mu, FILE *file)
{
    int iterations = 0;

    double x0 = left_end;
    double y0 = bottom_end;

    double x, y, nullity;

    while (1)
    {
        x = x0 - reverse_J_11(x0, y0) * f(x0, y0) - reverse_J_12(x0, y0) * g(x0, y0);
        y = y0 - reverse_J_21(x0, y0) * f(x0, y0) - reverse_J_22(x0, y0) * g(x0, y0);

        if (max(fabs(x - x0), fabs(y - y0)) < EPSILON / mu)
        {
            break;
        }

        x0 = x;
        y0 = y;
        iterations++;

        if (iterations > 10000)
        {
            printf ("\n Too many iterations !");
            fprintf(file, "\n Too many iterations !");
            return -2;
        }
    }

    nullity = max(fabs(f(x, y)), fabs(g(x, y)));

    printf ("\nNewton method gave the root (x, y) = (%f, %f) in %d iterations and ||f(x, y); g(x, y)|| = %f\n", x, y, iterations + 1, nullity);
    fprintf(file, "\nNewton method gave the root (x, y) = (%f, %f) in %d iterations and ||f(x, y); g(x, y)|| = %f\n", x, y, iterations + 1, nullity);
    return nullity;
}

//--------------------------------------------------------------------------------------------------------------------

int main(void)
{
    FILE* file;

    file = fopen("output.txt", "w");

    if (!file)
    {
        printf("Error opening file for writing. \n");
        return 0;
    }

    printf("Initial vector: (%f, %f)\n", left_end, bottom_end);
    fprintf(file, "Initial vector: (%f, %f)\n", left_end, bottom_end);

    double q  = J_phi_norm(file);
    double mu = mu_find   (file);

    simple_iteration_method(q, file);
    seidel_method          (file);

    printf ("\nmu = %f\n", mu);
    fprintf(file, "\nmu = %f\n", mu);

    newton_method          (mu, file);
    
    return 0;
}