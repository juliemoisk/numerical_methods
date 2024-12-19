#define _USE_MATH_DEFINES

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define alpha     3.0 // 1.0 2.0 3.0
#define EPSILON   1e-4
#define e         0.00005
#define left_end  1.0 //1.6 1.5 1.0
#define right_end 2.0


double function(double x) 
{
    return (1 + sin(x)) * sin(x) - alpha - 3 * x + 5;
}

double derivative(double x)
{
    return 2 * cos(x) * sin(x) + cos(x) - 3;
}

double second_derivative(double x)
{
    return -2*sin(x)*sin(x) - sin(x) + 2 * cos(x)*cos(x);
}

//--------------------------------------------------------------------------------------------------------------------

double phi(double x)
{
    return ((1 + sin(x)) * sin(x) - alpha + 5) / 3;
}

double phi_derivative(double x)
{
    return (2 * cos(x) * sin(x) + cos(x)) / 3;
}

double phi_second_derivative(double x)
{
    return (-2 * sin(x) * sin(x) - sin(x) + 2 * cos(x) * cos(x)) / 3;
}

//--------------------------------------------------------------------------------------------------------------------

int sign(double value)
{
    if (value > 0)
        return 1;
    else if (value < 0)
        return -1;
    else
        return 0;
}

void sign_on_the_ends(double (*func)(double), double a, double b, FILE *file)
{
    int sign_a = sign(func(a));
    int sign_b = sign(func(b));

    char signum_a[4];
    char signum_b[4];

    if      (sign_a > 0)  strcpy(signum_a, "+");
    else if (sign_a < 0)  strcpy(signum_a, "-");
    else                  strcpy(signum_a, "0");

    if      (sign_b > 0)  strcpy(signum_b, "+");
    else if (sign_b < 0)  strcpy(signum_b, "-");
    else                  strcpy(signum_b, "0");

    printf("\tsign for the left end:  %s\n", signum_a);
    printf("\tsign for the right end: %s\n", signum_b);

    fprintf(file, "\tsign for the left end:  %s\n", signum_a);
    fprintf(file, "\tsign for the right end: %s\n", signum_b);
}

//--------------------------------------------------------------------------------------------------------------------

double gradient_ascent(double (*func)(double), double (*derivative_f)(double), double a, double b, FILE *file)
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
            x1 += derivative_f(x1) * 0.01;
            temp = fabs(x1 - x2);

            if (function_max < func(x1) && x1 >= a && x1 <= b)
            {
                function_max = func(x1);
                x_max = x1;
            }
        }
    }
    printf("\tmaximum is in %.2f and equals %.2f\n", x_max, function_max);

    fprintf(file, "\tmaximum is in %.2f and equals %.2f\n", x_max, function_max);
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
            x1 -= derivative_f(x1) * 0.01;
            temp = fabs(x1 - x2);

            if (function_min > func(x1) && x1 >= a && x1 <= b) 
            {
                function_min = func(x1);
                x_min = x1;
            }
        }
    }

    printf("\tminimum is in %.2f and equals %.2f\n", x_min, function_min);

    fprintf(file, "\tminimum is in %.2f and equals %.2f\n", x_min, function_min);
    return function_min;
}

//--------------------------------------------------------------------------------------------------------------------

double chord_method(double m, double M, double a, double b, FILE *file) 
{
    int iterations = 0;
    double x0 = b;
    double x;

    if (function(a) * function(b) > 0) 
    {
        printf("ERROR ! Check the sign of the function in the ends of the interval ! \n");
        fprintf(file, "ERROR ! Check the sign of the function in the ends of the interval ! \n");
        return -1;
    }

    while (1)
    {
        x = x0 - function(x0) * (x0 - a) / (function(x0) - function(a));

        if (fabs(x - x0) < EPSILON * m / (M - m)) 
        {
            break;
        }

        x0 = x;
        iterations++;
        if (iterations > 10000) 
        {
            printf("Too many iterations\n");
            fprintf(file, "Too many iterations\n");
            return -2;
        }
    }

    printf("\nChord method gave the root x = %f in %d iterations and f(x) = %f\n", x, iterations + 1, function(x));
    fprintf(file, "\nChord method gave the root x = %f in %d iterations and f(x) = %f\n", x, iterations + 1, function(x));
    return x;
}

double newton_method(double m, double M, double a, double b, FILE *file)
{
    int iterations = 0;
    double x0 = b;
    double x;

    if (function(x0) * second_derivative(x0) <= 0)
    {
        printf("ERROR ! Not enough conditions for newton method convergence ! Check the 3rd condition");
        fprintf(file, "ERROR ! Not enough conditions for newton method convergence ! Check the 3rd condition");
        return -3;
    }

    if (derivative(a) * derivative(b) < 0)
    {
        printf("ERROR ! Function is not monotonous, not enough conditions for newton method convergence !");
        fprintf(file, "ERROR ! Function is not monotonous, not enough conditions for newton method convergence !");
        return -2;
    }

    if (second_derivative(a) * second_derivative(b) < 0)
    {
        printf("ERROR ! Derivative is not monotonous, not enough conditions for newton method convergence !");
        fprintf(file, "ERROR ! Derivative is not monotonous, not enough conditions for newton method convergence !");
        return -1;
    }

    while (1)
    {
        x = x0 - function(x0) / derivative(x0);

        if (fabs(x - x0) < EPSILON * m / (M - m))
        {
            break;
        }

        x0 = x;
        iterations++;
        if (iterations > 10000)
        {
            printf("Too many iterations\n");
            fprintf(file, "Too many iterations\n");
            return -2;
        }
    }

    printf("\nNewton method gave the root x = %f in %d iterations and f(x) = %f\n", x, iterations + 1, function(x));
    fprintf(file, "\nNewton method gave the root x = %f in %d iterations and f(x) = %f\n", x, iterations + 1, function(x));
    return x;
}

double simple_iteration_method(double M_phi, FILE *file)
{
    int iterations = 0;
    double x0 = right_end;
    double x;

    while (1)
    {
        x = phi(x0);

        if (fabs(x - x0) < EPSILON * (1 - M_phi) / M_phi)
        {
            break;
        }

        x0 = x;
        iterations++;
        if (iterations > 10000)
        {
            printf("Too many iterations\n");
            fprintf(file, "Too many iterations\n");
            return -2;
        }
    }

    printf("\nSimple iteration method gave the root x = %f in %d iterations and f(x) = %f\n", x, iterations + 1, function(x));
    fprintf(file, "\nSimple iteration method gave the root x = %f in %d iterations and f(x) = %f\n", x, iterations + 1, function(x));
    return x;
}

double bisection_method(double a, double b, FILE *file)
{
    double a1 = a, b1 = b;
    int n = int(log((b - a) / EPSILON) / log(2.0)) + 1; 
    double ksi = (a + b) / 2.0;
    int iterations = 0;

    for (int i = 0; i <= n && fabs(function(ksi)) > e; ++i)
    {
        if (function(ksi) * function(a1) < 0)
        {
            b1 = ksi;
        }
        else
        {
            a1 = ksi;
        }
        ksi = (a1 + b1) / 2.0;
        iterations++;
    }

    printf("\nBisection method gave the root ksi* = %f in %d iterations and f(ksi*) = %f\n", ksi, iterations, function(ksi));
    fprintf(file, "\nBisection method gave the root ksi* = %f in %d iterations and f(ksi*) = %f\n", ksi, iterations, function(ksi));
    return ksi;
}

double etkin_process(FILE *file)
{
    int iterations = 0;
    double x0 = right_end;
    double x1 = phi(x0);
    double x;

    while (1)
    {
        double x15 = phi(x1);
        x = (x0 * x15 - x1 * x1) / (x0 - 2 * x1 + x15);

        if (fabs(x - x1) < EPSILON)
        {
            break;
        }

        x0 = x1;
        x1 = x;

        iterations++;
        if (iterations > 10000)
        {
            printf("Too many iterations\n");
            fprintf(file, "Too many iterations\n");
            return -2;
        }
    }

    printf("\nEtkin process gave the root x = %f in %d iterations and f(x) = %f\n", x, iterations + 1, function(x));
    fprintf(file, "\nEtkin process gave the root x = %f in %d iterations and f(x) = %f\n", x, iterations + 1, function(x));
    return x;

}

//--------------------------------------------------------------------------------------------------------------------

int main(void)
{
    FILE *file;

    double derivative_min, derivative_max, m, M;
    double derivative_phi_min, derivative_phi_max, m_phi, M_phi;

    file = fopen("output.txt", "w");

    if (!file)
    {
        printf("Error opening file for writing. \n");
        return 0;
    }

    fprintf(file, "Analysis of function: f(x) = (1 + sin(x)) * sin(x) - alpha - 3 * x + 5\n\n");

    printf("Interval is [%.2f, %.2f], alpha = %.1f\n", left_end, right_end, alpha);
    fprintf(file,"Interval is [%.2f, %.2f], alpha = %.1f\n", left_end, right_end, alpha);
 
    printf("\nInitial function = (1 + sin(x)) * sin(x) - alpha - 3 * x + 5\n");
    sign_on_the_ends(function, left_end, right_end, file);

    printf("\nDerivative = 2 * cos(x) * sin(x) + cos(x) - 3\n");
    fprintf(file, "\nDerivative = 2 * cos(x) * sin(x) + cos(x) - 3\n");
    sign_on_the_ends(derivative, left_end, right_end, file);
    derivative_min = gradient_descent(derivative, second_derivative, left_end, right_end, file);
    derivative_max = gradient_ascent (derivative, second_derivative, left_end, right_end, file);

    printf("\nSecond derivative = -2*sin^2(x) - sin(x) + 2 * cos^2(x)\n");
    fprintf(file,"\nSecond derivative = -2*sin^2(x) - sin(x) + 2 * cos^2(x)\n");
    sign_on_the_ends(second_derivative, left_end, right_end, file);

    //---------------------------------------------------------------------------------------------------------------
    if (fabs(derivative_min) > fabs(derivative_max))
    {
        m = fabs(derivative_max);
        M = fabs(derivative_min);
    }
    else
    {
        m = fabs(derivative_min);
        M = fabs(derivative_max);
    }

    printf("\nm = %.2f\nM = %.2f\n", m, M);
    fprintf(file, "\nm = %.2f\nM = %.2f\n", m, M);

    //---------------------------------------------------------------------------------------------------------------

    printf("\nPhi = ((1 + sin(x)) * sin(x) - alpha + 5) / 3\n");
    fprintf(file,"\nPhi = ((1 + sin(x)) * sin(x) - alpha + 5) / 3\n");

    derivative_phi_min = gradient_descent(phi_derivative, phi_second_derivative, left_end, right_end, file);
    derivative_phi_max = gradient_ascent (phi_derivative, phi_second_derivative, left_end, right_end, file);

    if (fabs(derivative_phi_min) > fabs(derivative_phi_max))
    {
        m_phi = fabs(derivative_phi_max);
        M_phi = fabs(derivative_phi_min);
    }
    else
    {
        m_phi = fabs(derivative_phi_min);
        M_phi = fabs(derivative_phi_max);
    }

    printf("\nm_phi = %.2f\nM_phi = q = %.2f\n", m_phi, M_phi);
    fprintf(file, "\nm_phi = %.2f\nM_phi = q = %.2f\n", m_phi, M_phi);

    //---------------------------------------------------------------------------------------------------------------

    chord_method            (m, M, left_end, right_end, file);
    newton_method           (m, M, left_end, right_end, file);
    simple_iteration_method (M_phi, file);
    bisection_method        (left_end, right_end, file);
    etkin_process           (file);

    /*printf("\n Phi(%f) = %f, Phi(%f) = %f, Phi(pi/2) = %f", left_end, phi(left_end), right_end, phi(right_end), phi(M_PI / 2.0));*/ //для доказательства того, что фи на отрезке 1.6 - 2.0 есть сжим отобр
    fclose(file);
    return 0;
}

