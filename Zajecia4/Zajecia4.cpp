#include <iostream>
#include <functional>
#include "Funkcje.h"
using namespace std;

int iteracja = 200;
const int sizeofvector = 3;
double tolerancjabledu = 0.001;
double tolerancjareziduum = 0.001;

double* solveLinear(double* wektor)
{
	double * tmp = new double [3] {0.0, 0.0, 0.0};
	double Macierze[sizeofvector][sizeofvector + 1] =
	{
		{df1dx1(wektor[0]), df1dx2(wektor[1]), df1dx3(wektor[2]), fun1(wektor[0], wektor[1], wektor[2])},
		{df2dx1(wektor[0]), df2dx2(wektor[1]), df2dx3(wektor[2]), fun2(wektor[0], wektor[1], wektor[2])},
		{df3dx1(wektor[0]), df3dx2(wektor[1]), df3dx3(wektor[2]), fun3(wektor[0], wektor[1], wektor[2])}
	};

	for (int i = 0; i < sizeofvector; i++)
	{
		if (Macierze[i][i] == 0)
		{
			int k = i + 1;
			while (k < sizeofvector && Macierze[k][i] == 0)
			{
				k++;
			}
			if (k < sizeofvector)
			{
				for (int a = 0; a <= sizeofvector; a++)
				{
					std::swap(Macierze[i][a], Macierze[k][a]);
				}
			}
		}

		for (int j = i + 1; j < sizeofvector; j++)
		{
			double wspolczynnik = Macierze[j][i] / Macierze[i][i];
			for (int a = 0; a <= sizeofvector; a++)
			{
				Macierze[j][a] -= wspolczynnik * Macierze[i][a];
			}
		}
	}

	for (int i = sizeofvector - 1; i >= 0; i--)
	{
		double suma = 0;
		for (int j = i + 1; j < sizeofvector; j++)
		{
			suma += Macierze[i][j] * tmp[j];
		}
		tmp[i] = (Macierze[i][sizeofvector] - suma) / Macierze[i][i];
	}
	return tmp;
}


void odejmowanieWektorow(double* a, double* b)
{
	for (int i = 0; i < sizeofvector; i++)
	{
		a[i] -= b[i];
	}
}

double maxWektora(double* wektor)
{
	double result = wektor[0];
	for (int i = 1; i < sizeofvector; i++)
    {
        if(result < wektor[i])
            result = wektor[i];

    }
	return result;
}

void writeVector(double* wektor)
{
	for (int i = 0; i < sizeofvector; i++)
		cout << wektor[i] << endl;
	cout << endl;
}

void solveNonlinear(double* przyblizenie)
{
	double* delta = new double[3] {przyblizenie[0], przyblizenie[1], przyblizenie[2]};
	double poprzedni[3];
	double* estymator = new double[3] {0.0, 0.0, 0.0};
	double* reziduum = new double [3] {0.0, 0.0, 0.0};
	for (int i = 0; i < iteracja; i++)
	{
		delete[] delta;
		delta = solveLinear(przyblizenie);
		for (int i = 0; i < sizeofvector; i++)
			poprzedni[i] = przyblizenie[i];

		odejmowanieWektorow(przyblizenie, delta);

		for (int i = 0; i < sizeofvector; i++)
			estymator[i] = przyblizenie[i] - poprzedni[i];

		reziduum[0] = fun1(przyblizenie[0], przyblizenie[1], przyblizenie[2]);
		reziduum[1] = fun2(przyblizenie[0], przyblizenie[1], przyblizenie[2]);
		reziduum[2] = fun3(przyblizenie[0], przyblizenie[1], przyblizenie[2]);


		printf("Iteracja= %d, najwiekszy estymator= %lf, najwieksze reziduum= %lf\n", i, (abs(maxWektora(estymator))), (abs(maxWektora(reziduum))));
		cout << "Przyblizenie: \n";
		writeVector(przyblizenie);
		cout << "Delta: \n";
		writeVector(delta);
		cout << "------------------------" << endl;

		if ((abs(maxWektora(estymator))) <= tolerancjabledu && (abs(maxWektora(reziduum))) <= tolerancjareziduum)
		{
			cout << "Wynik z kryterium dokladnosci: "<< endl;
			writeVector(przyblizenie);
			return;
		}
	}
	cout << "Rozwiazano po " << iteracja << " iteracjach. Wynik to wektor:" << endl;
	writeVector(przyblizenie);
}

int main()
{
	double wektor[3] = { 10.0, 5.0, -3.0 };
	solveNonlinear(wektor);
	return 0;
}
