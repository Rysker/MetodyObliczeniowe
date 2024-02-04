#include <iostream>
#include <functional>
#include <cmath>
#include <Windows.h>
using namespace std;

int iteracje = 20;
double tolerancjabledu = 0.001;
double tolerancjareziduum = 0.001;

void Picard(function <double(double)> lambda, double x0)
{
    cout <<"Metoda Picarda"<< endl;
    int i = 1;
    double przyblizenie = x0;
    double estymator = 0;
    double reziduum = 0;
    double fix = 0;

    for(i; i <= iteracje; i++)
    {
        fix = przyblizenie + lambda(przyblizenie);
        estymator = fix - przyblizenie;
        przyblizenie = fix;
        reziduum = lambda(przyblizenie);
        printf("Iteracja %d, przyblizenie= %lf, estymator= %lf, reziduum= %lf\n", i, przyblizenie, estymator, reziduum);
        if(abs(estymator) <= tolerancjabledu && abs(reziduum) <= tolerancjareziduum)
        {
            cout <<"Wynik: " << przyblizenie << endl << endl;
            return;
        }
    }
    cout <<"Wynik z arbitralnym ograniczeniem iteracji: " << przyblizenie << endl << endl;
}

void bisekcja(function <double(double)> lambda, double a, double b)
{
    cout << "Metoda bisekcji" << endl;
    if((lambda(a) < 0 && lambda(b) < 0) || (lambda(a) >= 0 && lambda(b) >= 0))
    {
        cout << "Nie ma miejsc zerowych w podanym przedziale!" << endl;
        return;
    }
    //Funkcja musi byc ciagla i miec MZ w podanym przedziale
    double estymator = 0;
    double reziduum = 0;
    double srodek;
    int i = 1;
    double przyblizenie = (a + b)/2;

    for(i; i <= iteracje; i++)
    {
        srodek = (a + b)/2;
        if((lambda(a) < 0 && lambda(srodek) >= 0) || (lambda(a) >= 0 && lambda(srodek) < 0))
            b = srodek;
        if((lambda(srodek) < 0 && lambda(b) >= 0) || (lambda(srodek) >= 0 && lambda(b) < 0))
            a = srodek;
        przyblizenie = (a + b)/2;
        estymator = (b - a)/2;
        reziduum = lambda(przyblizenie);
        printf("Iteracja %d, przyblizenie= %lf, estymator= %lf, reziduum= %lf, a= %lf, b= %lf\n", i, przyblizenie, estymator, reziduum, a, b);

        if(abs(estymator) <= tolerancjabledu && abs(reziduum) <= tolerancjareziduum)
        {
            cout <<"Wynik: " << przyblizenie << endl << endl;
            return;
        }
    }
    cout <<"Wynik z arbitralnym ograniczeniem iteracji: " << przyblizenie << endl << endl;
}

void Newton(function <double(double)> lambda, function <double(double)> pochodna, double x0)
{
    cout << "Metoda Newtona" << endl;
    int i = 1;
    double tmp;
    double przyblizenie = x0;
    double estymator = 0;
    double reziduum = 0;

    for(i; i <= iteracje; i++)
    {
        tmp = przyblizenie;
        if(pochodna(przyblizenie) == 0)
        {
            cout << "Blad dzielenia przez 0! Nie mozna znalezc pierwiastka!" << endl;
            return;
        }
        przyblizenie = przyblizenie - (lambda(przyblizenie)/pochodna(przyblizenie));
        estymator = przyblizenie - tmp;
        reziduum = lambda(przyblizenie);
        printf("Iteracja %d, przyblizenie= %lf, estymator= %lf, reziduum= %lf\n", i, przyblizenie, estymator, reziduum);
        if(abs(estymator) <= tolerancjabledu && abs(reziduum) <= tolerancjareziduum)
        {
            cout <<"Wynik: " << przyblizenie << endl << endl;
            return;
        }
    }
    cout <<"Wynik z arbitralnym ograniczeniem iteracji: " << przyblizenie << endl << endl;
}

void sieczne(function <double(double)> lambda, double xn, double xn1)
{
    cout << "Metoda siecznych" << endl;
    if(xn == xn1)
    {
        cout << "Podano dwa te same punkty!" << endl;
        return;
    }

    int i = 1;
    double estymator = 0;
    double reziduum = 0;
    double przyblizenie = xn;

    for(i; i <= iteracje; i++)
    {
        double tmpxn1 = xn1;
        double tmpxn = xn;
        if((lambda(xn1)-lambda(xn)) == 0)
        {
            cout << "Blad dzielenia przez 0! Nie mozna znalezc pierwiastka!" << endl;
            return;
        }
        przyblizenie = xn1 - (lambda(xn1)*(xn1-xn))/(lambda(xn1)-lambda(xn));
        xn = xn1;
        xn1 = przyblizenie;
        estymator = xn1 - xn;
        reziduum = lambda(xn1);
        printf("Iteracja %d, przyblizenie= %lf, estymator= %lf, reziduum= %lf, xn1= %lf, xn= %lf\n", i, przyblizenie, estymator, reziduum, tmpxn1, tmpxn);
        if(abs(estymator) <= tolerancjabledu && abs(reziduum) <= tolerancjareziduum)
        {
            cout <<"Wynik: " << przyblizenie << endl << endl;
            return;
        }
    }
    cout <<"Wynik z arbitralnym ograniczeniem iteracji: " << przyblizenie << endl << endl;

}

int main()
{
    cout << "Funkcja 1" << endl;
    Picard([](double x)->double {return pow(sin(x/4), 2) - x;}, 3);
    bisekcja([](double x)->double {return pow(sin(x/4), 2) - x;}, -1, 1);
    Newton([](double x)->double {return pow(sin(x/4), 2) - x;}, [](double x)->double {return sin(x/2)/4 - 1;}, 2);
    sieczne([](double x)->double {return pow(sin(x/4), 2) - x;}, 5, 3);

    cout << "Funkcja 2" << endl;
    Picard([](double x)->double {return tan(2*x) - x - 1.0;}, 3);
    bisekcja([](double x)->double {return tan(2*x) - x - 1.0;}, -0.2, 0.5);
    Newton([](double x)->double {return tan(2*x) - x - 1.0;}, [](double x)->double {return 2/(pow(cos(2*x), 2));}, 0);
    sieczne([](double x)->double {return tan(2*x) - x - 1.0;}, 5, 3);
    return 0;
}
