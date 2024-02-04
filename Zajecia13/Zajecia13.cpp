#include <iostream>
#include <fstream>
#include "calerf.h"
using namespace std;

double erfc_value = 2.435813771487221;
double funkcja(double x)
{
    return calerfpack::erf_l(x);
}

double prostokaty_lewa(double x, long long unsigned podprzedzialy)
{
    double krok = (x - 0.) / podprzedzialy;
    double result = 0.;
    for (long long unsigned i = 0; i * krok < x; i++)
        result += funkcja(krok * i) * (krok);
    return result;
}

double prostokaty_prawa(double x, long long unsigned podprzedzialy)
{
    double krok = (x - 0.) / podprzedzialy;
    double result = 0.;
    for (long long unsigned i = 0; i * krok < x; i++)
        result += funkcja(krok * i + krok) * (krok);
    return result;
}

double prostokaty_srodek(double x, long long unsigned podprzedzialy)
{
    double krok = (x - 0.) / podprzedzialy;
    double result = 0.;
    for (long long unsigned i = 0; i * krok < x; i++)
    {
        double midpoint = (krok * (2 * i + 1)) / 2.;
        result += funkcja(midpoint) * (krok);
    }
    return result;
}

double trapezy(double x, long long unsigned podprzedzialy)
{
    double krok = (x - 0.) / podprzedzialy;
    double result = 0.;
    for (long long unsigned i = 0; i * krok < x; i++)
    {
        result += (funkcja(i * krok) + funkcja(i * krok + krok)) / 2. * (krok);
    }
    return result;
}

double parabole(double x, long long unsigned podprzedzialy)
{
    double krok = (x - 0.) / podprzedzialy;
    double result = 0.;
    for (long long unsigned i = 0; i * krok < x; i ++)
    {
        result += (1. / 6. * funkcja(i * krok) + 4. / 6. * funkcja(i * krok + krok / 2.) + 1. / 6. * funkcja(i * krok + krok)) * krok;
    }
    return result;
}


int main()
{
    ofstream fptr1, fptr2, fptr3, fptr4, fptr5;
    fptr1.open("prostokaty_lewa.txt");
    fptr2.open("prostokaty_prawa.txt");
    fptr3.open("prostokaty_srodek.txt");
    fptr4.open("trapezy.txt");
    fptr5.open("parabole.txt");
    double x = 3.;
    double przedzialy = 2.;
    double result;
    do
    {
        cout << "Przedzialy: " << przedzialy << endl;
        double krok = (x - 0.) / przedzialy;
        result = prostokaty_lewa(x, przedzialy);
        fptr1 << fabs((erfc_value - result)/erfc_value) << " " << krok << endl;

        result = prostokaty_prawa(x, przedzialy);
        fptr2 << fabs((erfc_value - result) / erfc_value) << " " << krok << endl;

        result = prostokaty_srodek(x, przedzialy);
        fptr3 << fabs((erfc_value - result) / erfc_value) << " " << krok << endl;

        result = trapezy(x, przedzialy);
        fptr4 << fabs((erfc_value - result) / erfc_value) << " " << krok << endl;

        result = parabole(x, przedzialy);
        fptr5 << fabs((erfc_value - result) / erfc_value) << " " << krok << endl;

        przedzialy *= 2;
    } while (przedzialy <= 1e8);
}