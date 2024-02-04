#include <iostream>
using namespace std;

const int sizeofvector = 4;
const int liczba_iteracji = 20;
double estymator_tolerancja[sizeofvector] = {0.1, 0.1, 0.1, 0.1};
double reziduum_tolerancja[sizeofvector] = {0.1, 0.1, 0.1, 0.1};
int indexes[sizeofvector];

void print_vector(double* v)
{
	for (int i = 0; i < sizeofvector; i++)
		cout << v[i] << " ";
}

void odejmij_wektory(double* v1, double* v2)
{
	for (int i = 0; i < sizeofvector; i++)
		v1[i] = v1[i] - v2[i];
}

bool check_estymator(double* v1, double* v2)
{
	double tmp[sizeofvector];
	for (int i = 0; i < sizeofvector; i++)
		tmp[i] = v1[i];

	odejmij_wektory(tmp, v2);

	for (int i = 0; i < sizeofvector; i++)
		tmp[i] = abs(tmp[i]);
	cout << "Estymator: ";
	print_vector(tmp);
	cout << endl;
	for (int i = 0; i < sizeofvector; i++)
	{
		if (tmp[i] > estymator_tolerancja[i])
			return false;
	}
	return true;
}

bool check_reziduum(double* v, double **macierz, double * b)
{
	double tmp[sizeofvector];
	
	for (int i = 0; i < sizeofvector; i++)
	{
		double suma = 0;
		for (int j = 0; j < sizeofvector; j++)
		{
			suma += macierz[i][j] * v[j];
		}
		tmp[i] = abs(suma - b[i]);
	}

	cout << "Reziduum: ";
	print_vector(tmp);
	cout << endl;
	for (int i = 0; i < sizeofvector; i++)
	{
		if (tmp[i] > reziduum_tolerancja[i])
			return false;
	}
	return true;
}

void Jacobi(double ** macierz, double* przyblizenie, double* b)
{
	cout << "Metoda Jacobiego" << endl;

	double u[sizeofvector] = { 0. };
	double xn[sizeofvector];
	for (int i = 0; i < sizeofvector; i++)
		xn[i] = przyblizenie[i];
	double tmp[sizeofvector];

	for (int k = 0; k < liczba_iteracji; k++)
	{
		cout << "Krok " << k << endl;
		for (int i = 0; i < sizeofvector; i++)
		{
			double suma = 0;
			for (int j = 0; j < sizeofvector; j++)
			{
				if (i != j)
					suma += macierz[i][j] * xn[j];
			}
			u[i] = (b[i] - suma) / macierz[i][i];
		}

		for (int i = 0; i < sizeofvector; i++)
		{
			tmp[i] = xn[i];
			xn[i] = u[i];
		}
		cout << "Przyblizenie: ";
		print_vector(xn);
		cout << endl;
		bool check1 = check_estymator(xn, tmp);
		bool check2 = check_reziduum(xn, macierz, b);
		if (check1 && check2)
		{
			cout << endl;
			break;
		}
		cout << "--------" << endl;
	}
}

void Gauss_Seidel(double** macierz, double* przyblizenie, double* b)
{
	cout << "Metoda Gaussa-Seidela" << endl;
	double u[sizeofvector] =  { 0 };
	for (int i = 0; i < sizeofvector; i++)
		u[i] = przyblizenie[i];
	double tmp[sizeofvector];
	for (int k = 0; k < liczba_iteracji; k++)
	{
		cout << "Krok " << k << endl;
		for (int i = 0; i < sizeofvector; i++)
		{
			double suma = 0;
			for (int j = 0; j < sizeofvector; j++)
			{
				if (i != j)
					suma += macierz[i][j] * u[j];
			}
			tmp[i] = u[i];
			u[i] = (b[i] - suma) / macierz[i][i];
		}

		cout << "Przyblizenie: ";
		print_vector(u);
		cout << endl;
		bool check1 = check_estymator(u, tmp);
		bool check2 = check_reziduum(u, macierz, b);
		if (check1 && check2)
		{
			cout << endl;
			break;
		}
		cout << "--------" << endl;
	}
}

void SOR(double** macierz, double* przyblizenie, double* b, double omega)
{
	cout << "Metoda SOR" << endl;
	double u[sizeofvector] = { 0 };
	for (int i = 0; i < sizeofvector; i++)
		u[i] = przyblizenie[i];
	double tmp[sizeofvector];

	for (int k = 0; k < liczba_iteracji; k++)
	{
		cout << "Krok " << k << endl;
		for (int i = 0; i < sizeofvector; i++)
		{
			double suma = 0;
			for (int j = 0; j < sizeofvector; j++)
			{
				if (i != j)
					suma += macierz[i][j] * u[j];
				else
					suma += (1. - 1. / omega) * macierz[i][i] * u[j];
			}
			tmp[i] = u[i];
			u[i] = (omega/macierz[i][i]) * (b[i] - suma);
		}

		cout << "Przyblizenie: ";
		print_vector(u);
		cout << endl;
		bool check1 = check_estymator(u, tmp);
		bool check2 = check_reziduum(u, macierz, b);
		if (check1 && check2)
		{
			cout << endl;
			break;
		}
		cout << "--------" << endl;
	}
}

int main()
{
	for (int i = 0; i < sizeofvector; i++)
		indexes[i] = i;

	double** A;
	A = new double* [4];
	A[0] = new double[4]{ 100., -1.0, 2.0, -3.0 };
	A[1] = new double[4]{ 1.0, 200.0, -4.0, 5.0 };
	A[2] = new double[4]{ -2.0, 4.0, 300.0, -6.0 };
	A[3] = new double[4]{ 3.0, -5.0, 6.0, 400.0 };

	double b[4] = { 116.0, -226.0, 912.0, -1174.0 };
	double przyblizenie[4] = { 2.0, 2.0, 2.0, 2.0 };

	Jacobi(A, przyblizenie, b);
	Gauss_Seidel(A, przyblizenie, b);
	SOR(A, przyblizenie, b, 0.5);
}