#include <iostream>
#include <cmath>
#include "calerf.h"
#include <fstream>

using namespace std;
//Wartości podane w treści zadania
double D = 1.;
double tmax = 2.;
double a = 6. * sqrt(D * tmax);
double x_pocz = -a;

//Warunki zakończenia dla metody iteracyjnej Gaussa-Seidela
int liczba_iteracji = 20;
double estymator_tolerancja = 1e-10;
double reziduum_tolerancja = 1e-10;

//Tablica z wartościami czasu, dla których przybliżone wyniki mają zostać zapisane do pliku
double przedzialy[5] = { 0.4, 0.9, 1.2, 1.5, 2.0 };

//Warunek poczatkowy równania różniczkowego
double warunek_poczatkowy(double x)
{
	if (x < 0)
		return 1.;
	else return 0.;
}

//Warunek brzegowy równania różniczkowego dla -inf
double warunek_brzegowy_1()
{
	return 1.;
}

//Warunek brzegowy równania różniczkowego dla +inf
double warunek_brzegowy_2()
{
	return 0.;
}

//Zwraca wartość rozwiązania analitycznego rozwiązywanego równania różniczkowego cząstkowego
double funkcja_analityczna(double x, double t)
{
	return (1./2.) * calerfpack::erfc_l(x/ (2.0 * sqrt(D * t)));
}

//Procedura obliczająca współczynniki ni na przekątnej
void thomasMatrix(double* u, double* d, double* l, int przekatna)
{
	for (int i = 1; i < przekatna; i++)
	{
		d[i] = d[i] - l[i - 1] * (1 / d[i - 1]) * u[i - 1];
	}
}

//Funkcja obliczająca ri w wektorze b i wyznaczająca rozwiązanie zgodnie z algorytmem Thomasa
double* thomasB(double* u, double* d, double* l, double* b, int przekatna)
{
	for (int i = 1; i < przekatna; i++)
	{
		b[i] = b[i] - l[i - 1] * (1 / d[i - 1]) * b[i - 1];
	}

	double* result = new double[przekatna];
	result[przekatna - 1] = (1 / d[przekatna - 1]) * b[przekatna - 1];
	for (int i = przekatna - 2; i >= 0; i--)
	{
		result[i] = (1 / d[i]) * (b[i] - (u[i] * result[i + 1]));
	}
	return result;
}

//Funkcja służąca do sprawdzenia estymatora dla otrzymanych x1, x2, ..., xn
//Zwraca true, jeśli każdy element z wektora v2 - v1 jest mniejszy niż, bądź równy tolerancji esytmatora
bool check_estymator(double* v1, double* v2, long long dl)
{
	double* tmp = new double[dl];

	tmp[0] = fabs(v1[0] - v2[0]);
	tmp[dl - 1] = fabs(v1[dl - 1] - v2[dl - 1]);

	for (int i = 1; i < dl - 1; i++)
	{
		tmp[i] = fabs(v1[i] - v2[i]);
	}

	for (int i = 0; i < dl; i++)
	{
		if (tmp[i] > estymator_tolerancja)
		{
			delete[] tmp;
			return false;
		}
	}

	delete[] tmp;
	return true;
}

//Funkcja służąca do sprawdzenia reziduum dla otrzymanego przybliżenia
//Zwraca true, jeśli każde wyliczane reziduum jest mniejsze bądź równe tolerancji reziduum
bool check_reziduum(double* v, double** macierz, double* b, long long dl)
{
	double* tmp = new double[dl];

	tmp[0] = fabs(macierz[0][1] * v[1] - (b[0] - macierz[0][0] * v[0]));
	tmp[dl - 1] = fabs(macierz[dl - 1][dl - 2] * v[dl - 2] - (b[dl - 1] - macierz[dl - 1][dl - 1] * v[dl - 1]));

	for (int i = 1; i < dl - 1; i++)
	{
		tmp[i] = fabs(macierz[i][i - 1] * v[i - 1] + macierz[i][i + 1] * v[i + 1] - (b[i] - macierz[i][i] * v[i]));
	}

	for (int i = 0; i < dl; i++)
	{
		if (tmp[i] > reziduum_tolerancja)
		{
			delete[] tmp;
			return false;
		}
	}

	delete[] tmp;
	return true;
}

//--------------------------------------------------------------------------------------------------------------
//Metoda iteracyjna Gaussa-Seidela służąca do rozwiązywania układów liniowych równań algebraicznych
//Wykorzystywana w metodzie Cranka-Nicolson do rozwiązywania powstających tam układów.
//--------------------------------------------------------------------------------------------------------------
void Gauss_Seidel(double** macierz, double* przyblizenie, double* b, long long dl)
{
	double* tmp = new double[dl];

	for (int k = 0; k < liczba_iteracji; k++)
	{
		for (int i = 0; i < dl; i++)
		{
			tmp[i] = przyblizenie[i];
			double suma = 0;

			if (i > 0)
				suma += macierz[i][i - 1] * przyblizenie[i - 1];

			if (i < dl - 1)
				suma += macierz[i][i + 1] * przyblizenie[i + 1];

			przyblizenie[i] = (b[i] - suma) / macierz[i][i];
		}

		bool check1 = check_estymator(przyblizenie, tmp, dl);
		bool check2 = check_reziduum(przyblizenie, macierz, b, dl);

		if (check1 && check2)
			break;
	}

	delete[] tmp;
}

//--------------------------------------------------------------------------------------------------------------
//Procedura Klasycznej Metody Bezpośredniej. Znajduje przybliżone rozwiązanie równania różniczkowego cząstkowego
//--------------------------------------------------------------------------------------------------------------
void klasyczna_metoda_bezposrednia()
{
	ofstream fptr1, fptr2;
	char nazwa_pliku[50];
	fptr1.open("bledy_krokt_KMB.txt");
	//Lambda podane w treści zadania. Współczynniki xh i th dobrane tak, aby spełniać równość lambda = th / (xh * xh)
	double lambda = 0.4;
	double xh = 0.005;
	double th = 0.00001;
	
	//Ilość punktów na jakie dzielimy nasz przedział [-a, a]
	unsigned long long dl = 2. * a / xh + 1;
	double* obecne, * poprzednie;
	obecne = new double[dl];
	poprzednie = new double[dl];

	//Wyliczenie wartości funkcji w czasie t = 0
	for (unsigned long long i = 0; i < dl; i++)
		poprzednie[i] = warunek_poczatkowy(i * xh - a);

	//W danej iteracji do while czas wynosi j * th
	long long j = 1;
	int zakres_czasowy = 0;
	
	//Obliczanie wartości funkcji w czasie t od th do 2
	do
	{
		//Wyliczenie wartości skrajnych wektora "obecne" z warunków brzegowych tj. x1, xn
		obecne[0] = warunek_brzegowy_1();
		obecne[dl - 1] = warunek_brzegowy_2();
		
		//Wyliczenie nowych wartości przybliżanego rozwiązania w punktach x2, x3, ..., xn-1
		for (unsigned long long i = 1; i < dl - 1; i++)
		{
			obecne[i] = lambda * poprzednie[i - 1] + (1. - 2. * lambda) * poprzednie[i] + lambda * poprzednie[i + 1];
		}
		
		//Zamieniamy wektory, tak aby obecne przybliżenie było poprzednim w następnej iteracji
		swap(poprzednie, obecne);

		//Wyliczanie największego błędu w danym czasie t, i zapisywanie go do pliku
		double blad = fabs(poprzednie[0] - funkcja_analityczna(- a, j * th));
		for (unsigned long long i = 1; i < dl; i++)
		{
			if (fabs(poprzednie[i] - funkcja_analityczna(i * xh - a, j * th)) > blad)
				blad = fabs(poprzednie[i] - funkcja_analityczna(i * xh - a, j * th));

		}
		fptr1 << j* th << " " << blad << endl;

		//Zapisywanie wartości uzyskanych przybliżeń, wraz z wartościami analitycznymi i błędami
		//Dla określonych wcześniej wartości czasu t
		if (j * th >= przedzialy[zakres_czasowy])
		{
			sprintf_s(nazwa_pliku, 50, "KMB%d.txt", zakres_czasowy);
			fptr2.open(nazwa_pliku);
			for (unsigned long long p = 0; p < dl; p++)
				fptr2 << p * xh - a << " " << poprzednie[p] << " " << funkcja_analityczna(p * xh - a, j * th) << " " << fabs(poprzednie[p] - funkcja_analityczna(p * xh - a, j * th)) << endl;
			fptr2.close();
			zakres_czasowy++;
		}
		//Zwiększenie czasu t
		j++;
			
	}while (j * th <= tmax);
}

//--------------------------------------------------------------------------------------------------------------
//Procedura Metody pośredniej Cranka-Nicolson.
//Parametr opcja:
//	-dla wartości 1 wykorzystuje algorytm Thomasa do rozwiązywania układów równań liniowych 
//	-dla wartości 2 wykorzystuje metodę iteracyjną Gaussa-Seidela do rozwiązywania układów równań liniowych 
//--------------------------------------------------------------------------------------------------------------
void Cranka_Nicolson(int opcja)
{
	ofstream fptr1, fptr2;
	char nazwa_pliku[50];

	//Metoda Cranka-Nicolson z użyciem algorytmu Thomasa
	if (opcja == 1)
	{
		//Plik przechowujący wartości błędów w funkcji czasu t
		fptr1.open("bledy_krokt_Thomas.txt");

		double lambda = 1.;
		double xh = 0.005;
		double th = 0.000025;

		//Ilość punktów na jakie dzielimy nasz przedział [-a, a]
		unsigned long long dl = 2. * a / xh + 1;
		double* d, * u, * l, * b, * poprzednie;
		d = new double[dl];
		u = new double[dl];
		l = new double[dl];
		b = new double[dl];
		poprzednie = new double[dl];
	
		//Wypełnienie wektora "poprzednie" za pomocą warunku początkowego
		for (long long i = 0; i < dl; i++)
			poprzednie[i] = warunek_poczatkowy(i * xh - a);
		
		//W danej iteracji do while czas wynosi i * th
		unsigned long long i = 1;
		int zakres_czasowy = 0;
		do
		{	//Wypełnienie "macierzy" A i wektora b układu równań liniowych Ax=b
			d[0] = 1.;
			u[0] = 0.;
			b[0] = 1;
			unsigned long long  j = 1;
			do
			{
				l[j - 1] = lambda / 2.;
				d[j] = -(1. + lambda);
				u[j] = lambda / 2.;
				b[j] = -(((lambda / 2.) * poprzednie[j - 1]) + ((1. - lambda) * poprzednie[j]) + (lambda / 2.) * (poprzednie[j + 1]));
				j++;
			}while (j < dl - 1);
			l[j - 1] = 0.;
			d[j] = 1.;
			b[j] = 0.;
			delete[] poprzednie;

			//Uzyskujemy rozwiązanie obecnego kroku za pomocą algorytmu Thomasa
			thomasMatrix(u, d, l, dl);
			poprzednie = thomasB(u, d, l, b, dl);

			//Znajdywanie największego błędu bezwzględnego w danym czasie t i zapisywanie do pliku
			double blad = fabs(poprzednie[0] - funkcja_analityczna(- a, i * th));
			for (long long q = 1; q < dl; q++)
			{
				if (fabs(poprzednie[q] - funkcja_analityczna(q * xh - a, i * th)) > blad)
				{
					blad = fabs(poprzednie[q] - funkcja_analityczna(q * xh - a, i * th));
				}
			}
			fptr1 << i*th << " " << blad << endl;

			//Zapisywanie wartości uzyskanych przybliżeń, wraz z wartościami analitycznymi i błędami
			//Dla określonych wcześniej wartości czasu t
			if (i * th >= przedzialy[zakres_czasowy])
			{
				sprintf_s(nazwa_pliku, 50, "Thomas%d.txt", zakres_czasowy);
				fptr2.open(nazwa_pliku);
				for (long long p = 0; p < dl; p++)
					fptr2 << p * xh - a << " " << poprzednie[p] << " " << funkcja_analityczna(p * xh - a, i * th) << " " << fabs(poprzednie[p] - funkcja_analityczna(p * xh - a, i * th)) << endl;
				fptr2.close();
				zakres_czasowy++;
			}
			//Zwiększamy czas
			i++;
		} while (i * th <= tmax);
			
		delete[] poprzednie;
		delete[] b;
		delete[] d;
		delete[] u;
		delete[] l;
	}

	//Metoda Cranka-Nicolson z użyciem metody iteracyjnej Gaussa-Seidela
	if (opcja == 2)
	{
		//Plik przechowujący wartości błędów w funkcji czasu t
		fptr1.open("bledy_krokt_Gauss.txt");
		double lambda = 1.;
		double xh = 0.005;
		double th = 0.000025;

		//Ilość punktów na jakie dzielimy nasz przedział [-a, a]
		unsigned long long dl = 2. * a / xh + 1;
		double * b, * poprzednie;
		double** macierz;

		//Tworzmy macierz, która będzie wykorzystywana w algorytmie Gaussa-Seidela
		macierz = new double* [dl];
		for (unsigned long long i = 0; i < dl; i++)
			macierz[i] = new double[dl] {0.};

		b = new double[dl];
		poprzednie = new double[dl];

		//Uzyskujemy rozwiązanie w czasie t=0 za pomocą warunku początkowego
		for (long long i = 0; i < dl; i++)
			poprzednie[i] = warunek_poczatkowy(i * xh - a);
		
		//W danej iteracji do while czas wynosi i * th
		unsigned long long i = 1;
		int zakres_czasowy = 0;
		do
		{
			//Wypełnienie macierzy A i wektora b układu równań liniowych Ax=b
			macierz[0][0] = 1.;
			b[0] = 1;
			unsigned long long  j = 1;
			do
			{
				macierz[j][j - 1] = lambda / 2.;
				macierz[j][j] = -(1. + lambda);
				macierz[j][j + 1] = lambda / 2.;
				b[j] = -(((lambda / 2.) * poprzednie[j - 1]) + ((1. - lambda) * poprzednie[j]) + (lambda / 2.) * (poprzednie[j + 1]));
				j++;
			}while (j < dl - 1);
			macierz[j][j] = 1.;
			b[j] = 0.;

			//Uzyskujemy rozwiązanie obecnego kroku za pomocą algorytmu Thomasa
			Gauss_Seidel(macierz, poprzednie, b, dl);

			//Znajdywanie największego błędu bezwzględnego w danym czasie t i zapisywanie do pliku
			double blad = fabs(poprzednie[0] - funkcja_analityczna(- a, i * th));
			for (unsigned long long q = 1; q < dl; q++)
			{
				if (fabs(poprzednie[q] - funkcja_analityczna(q * xh - a, i * th)) > blad)
					blad = fabs(poprzednie[q] - funkcja_analityczna(q * xh - a, i * th));

			}
			fptr1 << i * th << " " << blad << endl;


			//Zapisywanie wartości uzyskanych przybliżeń, wraz z wartościami analitycznymi i błędami
			//Dla określonych wcześniej wartości czasu t
			if (i * th >= przedzialy[zakres_czasowy])
			{
				sprintf_s(nazwa_pliku, 50, "Gauss%d.txt", zakres_czasowy);
				fptr2.open(nazwa_pliku);
				for (unsigned long long p = 0; p < dl; p++)
					fptr2 << p * xh - a << " " << poprzednie[p] << " " << funkcja_analityczna(p * xh - a, i * th) << " " << fabs(poprzednie[p] - funkcja_analityczna(p * xh - a, i * th)) << endl;
				fptr2.close();
				zakres_czasowy++;
			}
			//Zwiększamy czas
			i++;
		}while (i * th <= tmax);
			
		delete[] macierz;
		delete[] b;
		delete[] poprzednie;
	}
}

int main()
{
	klasyczna_metoda_bezposrednia();
	Cranka_Nicolson(1);
	Cranka_Nicolson(2);
}