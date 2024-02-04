#include <iostream>
using namespace std;

const int przekatna = 6;

void printArray(double* tab, int size)
{
	for (int i = 0; i < size; i++)
		cout << tab[i] << " ";
	cout << endl << "----------" << endl;
}

void thomasMatrix(double *u, double *d, double *l, int przekatna)
{
	for (int i = 1; i < przekatna; i++)
	{
		d[i] = d[i] - l[i - 1] * (1 / d[i - 1]) * u[i - 1];
	}
}

double* thomasB(double* u, double* d, double* l, double *b, int przekatna)
{
	for (int i = 1; i < przekatna; i++)
	{
		b[i] = b[i] - l[i - 1] * (1 / d[i - 1]) * b[i - 1];
	}

	cout << "Wektor r: ";
	printArray(b, przekatna);
	double* result = new double[przekatna];
	result[przekatna - 1] = (1 / d[przekatna - 1]) * b[przekatna - 1];
	for (int i = przekatna - 2; i >= 0; i--)
	{
		result[i] = (1 / d[i]) * (b[i] - (u[i] * result[i + 1]));
	}
	return result;
}

int main()
{
	double* u = new double[przekatna - 1]{ 1. / 2., 1. / 4., 1. / 6., 1. / 8., 1. / 10. };
	double* d = new double[przekatna] {10., 20., 30., 30., 20., 10.};
	double* l = new double[przekatna - 1]{ 1. / 3., 1. / 5., 1. / 7., 1. / 9., 1. / 11. };
	double* b = new double[przekatna] {31., 165. / 4., 917. / 30., 851. / 28., 3637. / 90., 332. / 11.};
	thomasMatrix(u, d, l, przekatna);

	cout << "Wektor ni: ";
	printArray(d, przekatna);

	double* ptr = thomasB(u, d, l, b, przekatna);
	for (int i = 0; i < przekatna; i++)
		cout << ptr[i] << endl;
}