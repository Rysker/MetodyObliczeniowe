#include <iostream>
using namespace std;

const int sizeofvector = 4;
int indexes[sizeofvector];

void decomposeA(double** Macierz)
{
	for (int i = 0; i < sizeofvector; i++)
	{
		if (Macierz[indexes[i]][i] == 0)
		{
			double max = 0;
			int index = i;
			int k = i + 1;
			for (int j = k; j < sizeofvector; j++)
			{
				double tmpmax = abs(Macierz[indexes[j]][i]);
				if (tmpmax != 0)
				{
					if (tmpmax > max)
					{
						max = tmpmax;
						index = j;
					}
				}
			}

			std::swap(indexes[i], indexes[index]);
		}

		for (int j = i + 1; j < sizeofvector; j++)
		{
			double wspolczynnik = Macierz[indexes[j]][i] / Macierz[indexes[i]][i];

			for (int a = i; a <= sizeofvector; a++)
			{
				Macierz[indexes[j]][a] -= wspolczynnik * Macierz[indexes[i]][a];
			}
			Macierz[indexes[j]][i] = wspolczynnik;
		}

		for (int i = 0; i < sizeofvector; i++)
		{
			for (int j = 0; j < sizeofvector; j++)
				printf("%.2f\t", Macierz[indexes[i]][j]);
			cout << endl;
		}
		cout << "----------------------------" << endl;
	}
}


double* solveLinear(double** Macierz, double* b)
{
	double y[sizeofvector] = { 0 };
	double* x = new double[sizeofvector] { 0 };
	// solve Ly = b
	for (int i = 0; i < sizeofvector; i++) {
		y[i] = b[indexes[i]];
		for (int j = 0; j < i; j++) {
			y[i] -= Macierz[indexes[i]][j] * y[j];
		}
	}

	// solve Ux = y
	for (int i = sizeofvector - 1; i >= 0; i--) {
		x[i] = y[i];
		for (int j = i + 1; j < sizeofvector; j++) {
			x[i] -= Macierz[indexes[i]][j] * x[j];
		}
		x[i] /= Macierz[indexes[i]][i];
	}

	return x;
}


int main()
{
	for (int i = 0; i < sizeofvector; i++)
		indexes[i] = i;

	double** A;
	A = new double* [4];
	A[0] = new double[4]{ 1.0, -20.0, 30.0, -4.0 };
	A[1] = new double[4]{ 2.0, -40.0, -6.0, 50.0 };
	A[2] = new double[4]{ 9.0, -180.0, 11.0, -12.0 };
	A[3] = new double[4]{ -16.0, 15.0, -140.0, 13.0 };

	double b[4] = { 35.0, 104.0, -366.0, -354.0 };

	decomposeA(A);
	for (int i = 0; i < sizeofvector; i++)
	{
		for (int j = 0; j < sizeofvector; j++)
			printf("%.2f\t", A[i][j]);
		cout << endl;
	}
	cout << endl << "Indexes: ";
	for (int j = 0; j < sizeofvector; j++)
	{
		printf("%d\t", indexes[j]);
	}
	cout << endl << "Wektor: ";
	double* ptr = solveLinear(A, b);
	for (int i = 0; i < sizeofvector; i++)
		printf("%.2f\t", ptr[indexes[i]]);
	cout << endl << "\tRozwiazanie\t" << endl;
	for (int i = 0; i < sizeofvector; i++)
		printf("%.2f\t", ptr[i]);
	cout << endl;
}
