#define _USE_MATH_DEFINES

#include <iostream>
#include <cmath>
#include <iostream>
#include <fstream>

using namespace std;
const int iteracje = 20;
double alfa = 0., beta = 1., gamma = -1., fi = 0., psi = 1., theta = 0.;
double x_konc = 100.;

double funkcja(double t, double y)
{
    return -((10. * t * t + 20.) / (t * t + 1.)) * (y -1.);
}

double funkcja2(double t)
{
    return 1. - exp(-10. *(t + atan(t)));
}

double funkcja3(double y_obecny, double y_poprzedni, double delta_t, double t)
{
    return (y_obecny - y_poprzedni) / delta_t - funkcja(t, y_obecny);
}

double funkcja4(double y_obecny, double y_poprzedni, double delta_t, double t)
{
    return (y_obecny - y_poprzedni) / delta_t - (funkcja(t - delta_t, y_poprzedni) + funkcja(t, y_obecny))/2.;
}

void Euler()
{
    ofstream fptr1, fptr2;
    //fptr1.open("wynik1.txt");
    fptr2.open("wynik4f.txt");
    double delta_t = 0.2;
    while (delta_t > 1e-4)
    {
        cout << "delta t= " << delta_t << endl;
        int dl = x_konc / delta_t + 1;
        double y_poprzedni = 0.;
        double y_obecny = 0.;
        double maxblad = fabs(y_obecny - funkcja2(0.));
        double i = delta_t;
        double j = 1;
        do
        {
            fptr2 << i - delta_t << " " << y_poprzedni << endl;
            y_obecny = y_poprzedni + funkcja(i, y_poprzedni) * delta_t;
            if (fabs(y_obecny - funkcja2(j * delta_t)) > maxblad)
                maxblad = fabs(y_obecny - funkcja2(j * delta_t));

            y_poprzedni = y_obecny;
            j++;
            i += delta_t;
        } while (j * delta_t < x_konc);

        //fptr1 << log10(delta_t) << " " << log10(maxblad) << endl;
        delta_t /= 2.;
    }
    //fptr1.close();
}

double sieczne(double y_poprzedni, double delta_t, double t, int opcja)
{
    int xn1 = -1;
    int xn = 0.5;
    double tolerancjabledu = 1e-4;
    double tolerancjareziduum = 1e-4;
    //cout << "Metoda siecznych" << endl;
    if (xn == xn1)
        return 0.;
    int i = 1;
    double estymator = 0;
    double reziduum = 0;
    if (opcja == 1)
    {
        double przyblizenie = y_poprzedni + funkcja(i, y_poprzedni) * delta_t;

        for (i; i <= iteracje; i++)
        {
            double tmpxn1 = xn1;
            double tmpxn = xn;
            przyblizenie = (xn1 - funkcja3(xn1, y_poprzedni, delta_t, t) * (xn1 - xn)) / (funkcja3(xn1, y_poprzedni, delta_t, t) - funkcja3(xn, y_poprzedni, delta_t, t));
            xn = xn1;
            xn1 = przyblizenie;
            estymator = xn1 - xn;
            reziduum = funkcja3(przyblizenie, y_poprzedni, delta_t, t);
            //printf("Iteracja %d, przyblizenie= %lf, estymator= %lf, reziduum= %lf, xn1= %lf, xn= %lf\n", i, przyblizenie, estymator, reziduum, tmpxn1, tmpxn);
            if (abs(estymator) <= tolerancjabledu && abs(reziduum) <= tolerancjareziduum)
                return przyblizenie;
        }
        return przyblizenie;
    }

    if (opcja == 2)
    {
        double przyblizenie = y_poprzedni + funkcja(i, y_poprzedni) * delta_t;

        for (i; i <= iteracje; i++)
        {
            double tmpxn1 = xn1;
            double tmpxn = xn;
            if (funkcja4(xn1, y_poprzedni, delta_t, t) - funkcja4(xn, y_poprzedni, delta_t, t) == 0)
                return 0.;
            przyblizenie = (xn1 - funkcja4(xn1, y_poprzedni, delta_t, t) * (xn1 - xn)) / (funkcja4(xn1, y_poprzedni, delta_t, t) - funkcja4(xn, y_poprzedni, delta_t, t));
            xn = xn1;
            xn1 = przyblizenie;
            estymator = xn1 - xn;
            reziduum = funkcja4(przyblizenie, y_poprzedni, delta_t, t);
            //printf("Iteracja %d, przyblizenie= %lf, estymator= %lf, reziduum= %lf, xn1= %lf, xn= %lf\n", i, przyblizenie, estymator, reziduum, tmpxn1, tmpxn);
            if (abs(estymator) <= tolerancjabledu && abs(reziduum) <= tolerancjareziduum)
                return przyblizenie;
        }
        return przyblizenie;
    }
}

void PosredniaEuler()
{
    ofstream fptr1, fptr2;
    //fptr1.open("wynik2.txt");
    fptr2.open("wynik2f.txt");
    double delta_t = 1e-5;
    while (delta_t > 3e-6)
    {
        cout << "delta t= " << delta_t << endl;
        double y_poprzedni = 0.;
        double y_obecny = 0.;
        double maxblad = fabs(y_obecny - funkcja2(0.));
        double i = delta_t;
        double j = 1;
        do
        {
            fptr2 << i - delta_t << " " << y_poprzedni << endl;
            y_obecny = sieczne(y_poprzedni, delta_t, i, 1);
            if (fabs(y_obecny - funkcja2(i)) > maxblad)
                maxblad = fabs(y_obecny - funkcja2(i));
            y_poprzedni = y_obecny;
            j++;
            i += delta_t;
        } while (i < x_konc);

        //fptr1 << log10(delta_t) << " " << log10(maxblad) << endl;
        delta_t /= 2.;
    }
    //fptr1.close();
}

void PosredniaTrapezow()
{
    ofstream fptr1, fptr2;
    //fptr1.open("wynik3.txt");
    fptr2.open("wynik3f.txt");
    double delta_t = 1e-5;
    while (delta_t > 3e-6)
    {
        cout << "delta t= " << delta_t << endl;
        double y_poprzedni = 0.;
        double y_obecny = 0.;
        double maxblad = fabs(y_obecny - funkcja2(0.));
        double i = delta_t;
        double j = 1;
        do
        {
            fptr2 << i - delta_t << " " << y_poprzedni << endl;
            y_obecny = sieczne(y_poprzedni, delta_t, i, 2);
            if (fabs(y_obecny - funkcja2(i)) > maxblad)
            {
                maxblad = fabs(y_obecny - funkcja2(i));
            }
            y_poprzedni = y_obecny;
            j++;
            i += delta_t;
        } while (i < x_konc);

        //fptr1 << log10(delta_t) << " " << log10(maxblad) << endl;
        delta_t /= 2.;
    }
    //fptr1.close();
}

int main()
{
    Euler();
   //PosredniaEuler();
   //PosredniaTrapezow();
}
