#ifndef INC_4_MINIMISATION_H
#define INC_4_MINIMISATION_H

#include <vector>
#include <iostream>
#include <cmath>

//calcule la fonction à l'étude
double f(std::vector<double> x);

//trouve le n minimum pour une précison donnée (pour Fibonacci)
void recherche(const std::vector<double> &a, std::vector<double> &b, double eps, long long int &fib1,
               long long int &fib2,
               long long int &fib3, long long int &n);

//fonctions pour les vecteurs
std::vector<double> operator+(const std::vector<double> &v1, const std::vector<double> &v2);

std::vector<double> operator-(const std::vector<double> &v1, const std::vector<double> &v2);

std::vector<double> operator*(const std::vector<double> &v, double nombre);

void demontrer(const std::vector<double> &v);

double distance(const std::vector<double> &v1, const std::vector<double> &v2);

std::vector<double> gradient(const std::vector<double> &x);

//trouve le minimum de la fonction multidimensionelle
std::vector<double> minMultidim(std::vector<double> x, double eps1, double eps2);

//trouve l'intervalle du minimum (sur une droite)
std::vector<double> minimisation(std::vector<double> x, const std::vector<double> &s, double eps);

//trouve le minimum sur l'intervalle
std::vector<double> fibonacci(std::vector<double> &a, std::vector<double> &b, double eps);


std::vector<double> gradient(const std::vector<double> &x)
{
    std::vector<double> s(x.size());
    s[0] = 600*x[0]*x[0]*x[0]*x[0]*x[0]-600*x[0]*x[0]*x[1]+2*x[0]-2; //-600 * x[0]*x[0] *(- x[0]*x[0]*x[0] + x[1]) + 2*x[0]-2;
    s[1] = -200 * x[0] * x[0] * x[0] + 200 * x[1];
    std::cout <<std::endl<< "Gradient. "<<std::endl;
    demontrer(x);
    demontrer(s);
    return s;
}

std::vector<double> minMultidim(std::vector<double> x, double eps1, double eps2)
{
    std::vector<double> xPrec;
    bool term = false;
    int pas = 0;
    while (!term)
    {
        ++pas;
        xPrec = x;
        x = minimisation(x, gradient(x), eps1 * eps2);
        //std::cout << "coucou " << f(x) << " " << f(xPrec) << " " << xPrec[0] << " " << xPrec[1] << " " << x[0] << " "
        //          << x[1] << std::endl << std::endl;
        std::cout<<"!"<<f(xPrec) - f(x)<<std::endl;
        if (fabs(f(xPrec) - f(x)) < eps1)
            break;
        term = true;
        for (long long int i = 0; i < x.size(); ++i)
            if (fabs(xPrec[i] - x[i]) > eps2)
            {
                term = false;
                break;
            }
    }
    std::cout <<std::endl<< "Fin réussie. " << "Nombre de pas: " << pas << ". Valeur : "<<f(x)<<std::endl;
    return x;
}

void
recherche(const std::vector<double> &a, std::vector<double> &b, double eps, long long int &fib1, long long int &fib2,
          long long int &fib3, long long int &n)
{
    double compar = distance(a, b) / eps;

    fib1 = 1, fib2 = 1, fib3 = 2, n = 1;
    while ((fib3) < compar)
    {
        ++n;
        fib1 = fib2;
        fib2 = fib3;
        fib3 = fib2 + fib1;
    }
}

std::vector<double> fibonacci(std::vector<double> &a, std::vector<double> &b, double eps)
{
    long long int n, fib1, fib2, fib3, pas = 0;
    recherche(a, b, eps, fib1, fib2, fib3, n);
    std::vector<double> x1, x2;
    double y1, y2;

    while (n > 1)
    {
        ++pas;
        x1 = a + (b - a) * (fib1 / (double) fib3);
        x2 = a + (b - a) * (fib2 / (double) fib3);
        y1 = f(x1);
        y2 = f(x2);
        if (y1 > y2)
            a = x1;
        else
            b = x2;
        --n;
        fib3 = fib2;
        fib2 = fib1;
        fib1 = fib3 - fib2;
        //demontrer(x1);
        //demontrer(x2);
    }
    std::cout << "Fibonacci. " << "Nombre de pas: " << pas << std::endl;
    demontrer((x1 + x2) * (0.5));
    return (x1 + x2) * (0.5);
}


double f(std::vector<double> x)
{
    return (1 - x[0]) * (1 - x[0]) + 100 * (x[1] - x[0]*x[0]*x[0]) * (x[1] - x[0]*x[0]*x[0]);
}


std::vector<double> operator+(const std::vector<double> &v1, const std::vector<double> &v2)
{
    std::vector<double> rep(v1.size());
    for (long long int i = 0; i < v1.size(); ++i)
        rep[i] = v1[i] + v2[i];
    return rep;
}

std::vector<double> operator-(const std::vector<double> &v1, const std::vector<double> &v2)
{
    std::vector<double> rep(v1.size());
    for (long long int i = 0; i < v1.size(); ++i)
        rep[i] = v1[i] - v2[i];
    return rep;
}

std::vector<double> operator*(const std::vector<double> &v, double nombre)
{
    std::vector<double> rep(v.size());
    for (long long int i = 0; i < v.size(); ++i)
        rep[i] = v[i] * nombre;
    return rep;
}

void demontrer(const std::vector<double> &v)
{
    for (const auto &x:v)
    {
        std::cout << x << ' ';
    }
    std::cout << std::endl;
}

double distance(const std::vector<double> &v1, const std::vector<double> &v2)
{
    double res = 0;
    for (long long int i = 0; i < v1.size(); ++i)
        res += (v1[i] - v2[i]) * (v1[i] - v2[i]);
    return sqrt(res);
}


std::vector<double> minimisation(std::vector<double> x, const std::vector<double> &s, double eps)
{
    long long int pas = 1;
    std::vector<double> prec(x.size()), suiv(x.size());
    double h = 1;
    bool term = false;

    prec = x - s;
    suiv = x + s;
    //std::cout<<f(x-s*1000) <<" "<< f(prec)<<" ";
    //demontrer(x-s*1000);
    demontrer(prec);
    if (f(x) > f(prec))
        if (f(x) > f(suiv))
            std::cout << "Err";
        else
        {
            h *= -1;
            swap(prec, suiv);
        }

    else if (f(x) < f(suiv))
        term = true;

    if (!term)
        do
        {
            ++pas;
            h *= 2;
            prec = x;
            x = suiv;
            suiv = suiv + s * h;
        } while (f(x) > f(suiv));

    std::cout << "Intervalle trouvée. " << "Nombre de pas: " << pas << std::endl;
    demontrer(prec);

    demontrer(suiv);
    return (fibonacci(prec, suiv, eps));
}

#endif //INC_4_MINIMISATION_H
