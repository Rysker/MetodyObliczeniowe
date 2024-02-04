#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

using namespace std;
const double PI = 3.1415926535897932384626433832795028841971693993751058209;

template<typename T>
T roznica_progresywna(T x, T h)
{
    return (sin(x + h) - sin(x)) / h;
}

template<typename T>
T roznica_wsteczna(T x, T h)
{
    return (sin(x) - sin(x - h)) / h;
}

template<typename T>
T roznica_centralna(T x, T h)
{
    return (sin(x + h) - sin(x - h)) / (2 * h);
}

template<typename T>
T roznica_3punktowapoczatek(T x, T h)
{
    T a = static_cast<T>(-3. / 2.);
    T b = static_cast<T>(2.);
    T c = static_cast<T>(-1. / 2.);
    return (a * sin(x) + b * sin(x + h) + c * sin(x + 2 * h)) / h;
}

template<typename T>
T roznica_3punktowakoniec(T x, T h)
{
    T a = static_cast<T>(1. / 2.);
    T b = static_cast<T>(2.);
    T c = static_cast<T>(3. / 2.);
    return (a * sin(x - 2 * h) - 2 * sin(x - h) + c * sin(x)) / h;
}

template<typename T>
void pochodne(T x1, T x2, T x3)
{
    ofstream fptr1, fptr2, fptr3, fptr4, fptr5, fptr6, fptr7;
    T f1, f2, f3;
    fptr1.open("wynik1.txt");
    fptr2.open("wynik2.txt");
    fptr3.open("wynik3.txt");
    fptr4.open("wynik4.txt");
    fptr5.open("wynik5.txt");
    fptr6.open("wynik6.txt");
    fptr7.open("wynik7.txt");
    T tab[3] = { static_cast<T>(0.5), static_cast<T>(0.2), static_cast<T>(0.1) };
    do
    {
        T h = tab[0];
        for (int i = 0; i < 3; i++)
        {
            if (h == 0.5)
            {
                h = tab[2];
                i = 2;
            }
            else
                h = tab[i];
            cout << "h= " << h << endl;
            cout << "x1 = " << x1 << endl;
            f1 = roznica_progresywna(x1, h);
            f2 = roznica_3punktowapoczatek(x1, h);
            cout << "Przyblizenie dwupunktowe: " << f1 << " " << ", trzypunktowe: " << f2 << endl << endl;
            fptr1 << log10(h) << " " << log10(fabs(f1 - cos(x1))) << endl;
            fptr2 << log10(h) << " " << log10(fabs(f2 - cos(x1))) << endl;
            cout << "x2 = " << x2 << endl;
            f1 = roznica_progresywna(x2, h);
            f2 = roznica_wsteczna(x2, h);
            f3 = roznica_centralna(x2, h);
            cout << "Przyblizenie dwupunktowe (roznica progresywna): " << f1 << " " << ", przyblizenie dwupunktowe (roznica wsteczna): " << f2 << " " << ", przyblizenie dwupunktowe (roznica centralna)/trzypunktowe: " << f3 << endl << endl;
            fptr3 << log10(h) << " " << log10(fabs(f1 - cos(x2))) << endl;
            fptr4 << log10(h) << " " << log10(fabs(f2 - cos(x2))) << endl;
            fptr5 << log10(h) << " " << log10(fabs(f3 - cos(x2))) << endl;
            cout << "x3 = " << x3 << endl;
            f1 = roznica_wsteczna(x3, h);
            f2 = roznica_3punktowakoniec(x3, h);
            cout << "Przyblizenie dwupunktowe: " << f1 << " " << ", trzypunktowe: " << f2 << endl << endl;
            fptr6 << log10(h) << " " << log10(fabs(f1 - cos(x3))) << endl;
            fptr7 << log10(h) << " " << log10(fabs(f2 - cos(x3))) << endl;
            cout<<"----------------------------------------------------------" << endl;
        }
        for (int i = 0; i < 3; i++)
            tab[i] = tab[i] / 10;
        h = tab[0];
    } while (tab[2] > 10e-18);
}

int main()
{
    float x1 = 0.f;
    float x2 = float(PI / 4.f);
    float x3 = float(PI / 2.f);
    pochodne(x1, x2, x3);

    double x4 = 0.;
    double x5 = PI / 4.;
    double x6 = PI / 2.;
    //pochodne(x4, x5, x6);
}
