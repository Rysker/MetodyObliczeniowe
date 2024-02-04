#include <iostream>
#include <fstream>
#include <cmath>
#define M_PI  3.14159265358979323846  

using namespace std;
double x_pocz = -1.;
double x_konc = 1.;

double funkcja(double x)
{
	return 1. / (1. + (10. * x * x * x * x * x * x));
}

void bazaNewton(double* x, double* y, int liczba_wezlow, char typ_wezlow)
{
	ofstream fptr1;
	char nazwa_pliku[50];
	sprintf_s(nazwa_pliku, 50, "%c%d.txt", typ_wezlow, liczba_wezlow);
	fptr1.open(nazwa_pliku);
	int n_liczba = 400;
	double** c;

	c = new double* [liczba_wezlow];

	for (int i = 0; i < liczba_wezlow; i++)
	{
		c[i] = new double[liczba_wezlow] {0.};
		c[i][0] = y[i];
	}

	for (int i = 1; i < liczba_wezlow; ++i)
	{
		for (int j = 0; j < liczba_wezlow - i; ++j)
		{
			c[j][i] = (c[j + 1][i - 1] - c[j][i - 1]) / (x[i + j] - x[j]);
		}
	}

	double krok2 = (x_konc - x_pocz) / (n_liczba - 1.);
	for (int i = 0; i < n_liczba; i++)
	{
		double wynik = c[0][liczba_wezlow - 1];
		for (int j = liczba_wezlow - 1; j > 0; --j) 
		{
			wynik *= ((i * krok2 + x_pocz) - x[j - 1]);
			wynik += c[0][j - 1];
		}
		fptr1 << (i * krok2 + x_pocz) << " " << wynik << endl;
	}
	fptr1.close();
}

void Czebyszew(int liczba_wezlow) 
{
	double* x, * y;
	double krok = (x_konc - x_pocz) / (liczba_wezlow - 1.);
	x = new double[liczba_wezlow];
	y = new double[liczba_wezlow];
	for (int i = 0; i < liczba_wezlow; i++) 
	{
		double pierwiastek = cos(((2. * i + 1.) / (2.0 * (liczba_wezlow - 1.0) + 2.0)) * M_PI);
		x[i] = (x_konc + x_pocz) / 2. + (x_konc - x_pocz) / 2. * pierwiastek;
		y[i] = funkcja(x[i]);
	}
	bazaNewton(x, y, liczba_wezlow, 'c');
}

void rownoodlegle(int liczba_wezlow)
{
	double* x, * y;
	double krok = (x_konc - x_pocz) / (liczba_wezlow - 1.);
	x = new double[liczba_wezlow];
	y = new double[liczba_wezlow];
	for (int i = 0; i < liczba_wezlow; i++)
	{
		x[i] = i * krok + x_pocz;
		y[i] = funkcja(x[i]);
	}
	bazaNewton(x, y, liczba_wezlow, 'r');
}

int main() {
	int wezly1 = 4;
	int wezly2 = 6;
	int wezly3 = 10;
	int wezly4 = 15;
	rownoodlegle(wezly1);
	rownoodlegle(wezly2);
	rownoodlegle(wezly3);
	rownoodlegle(wezly4);
	Czebyszew(wezly1);
	Czebyszew(wezly2);
	Czebyszew(wezly3);
	Czebyszew(wezly4);
}
