#define _USE_MATH_DEFINES

#include <iostream>
#include <cmath>
#include <iostream>
#include <fstream>

using namespace std;

double alfa = 0., beta = 1., gamma = -1., fi = 0., psi = 1., theta = 0.;

double funkcja(double x)
{
	return (pow(M_E, 2 - 2 * x) - 4 * pow(M_E, 4 - 2 * x) + 4 * pow(M_E, 2 * x) - pow(M_E, 2 + 2 * x) - x + x * pow(M_E, 4)) / (4 - 4 * pow(M_E, 4));
}

void printArray(double* tab, int size)
{
	for (int i = 0; i < size; i++)
		cout << tab[i] << " ";
	cout << endl << "----------" << endl;
}

void thomasMatrix(double* u, double* d, double* l, int przekatna)
{
	for (int i = 1; i < przekatna; i++)
	{
		d[i] = d[i] - l[i - 1] * (1 / d[i - 1]) * u[i - 1];
	}
}

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

void konwencjonalna()
{
    ofstream fptr1, fptr2;
    //fptr1.open("wynik1.txt");
    fptr2.open("wynik1f.txt");
    double h = 1e-4;
    while (h > 1e-5)
    {
        cout << "h= " << h << endl;
        int dl = 1.0 / h + 1;
        double* d, * u, * l, * b;
        d = new double[dl];
        u = new double[dl];
        l = new double[dl];
        b = new double[dl];
        d[0] = (beta - alfa / h);
        u[0] = alfa / h;
        b[0] = -gamma;
        int j = 1;
        double i = h;
        do
        {
            l[j - 1] = 1. / (h * h);
            d[j] = -4. - 2. / (h * h);
            u[j] = 1. / (h * h);
            b[j] = i;
            j++;
            i += h;
        } while (j < dl - 1);

        l[j - 1] = -fi / h;
        d[j] = (fi / h + psi);
        b[j] = -theta;
        thomasMatrix(u, d, l, dl);

        double* ptr = thomasB(u, d, l, b, dl);
        double maxblad = fabs(ptr[0] - funkcja(0));
        for (int i = 0; i < dl; i++)
        {
            if (fabs(ptr[i] - funkcja(0. + i * h)) > maxblad)
                maxblad = fabs(ptr[i] - funkcja(0. + i * h));
        }
        for (int i = 0; i < dl; i++)
            fptr2 << i * h << " " << ptr[i] << endl;
        //fptr1 << log10(h) << " " << log10(maxblad) << endl;
        delete[] d;
        delete[] u;
        delete[] l;
        delete[] b;
        delete[] ptr;
        h /= 2.;
    }
    //fptr1.close();
}

void Numerow()
{
    ofstream fptr1, fptr2;
    fptr1.open("wynik2tmp.txt");
    //fptr2.open("wynik2f.txt");
    double h = 0.01;
    //double h = 1e-4;
    while (h > 1e-8)
    {
        cout << "h= " << h << endl;
        int dl = 1.0 / h + 1;
        double* d, * u, * l, * b;
        d = new double[dl];
        u = new double[dl];
        l = new double[dl];
        b = new double[dl];
        d[0] = beta - alfa / h;
        u[0] = alfa / h;
        b[0] = -gamma;
        int j = 1;
        double i = h;
        do
        {
            l[j - 1] = 1. / (h * h) - 1. / 3.;
            d[j] = - 10. / 3. - 2. /(h*h);
            u[j] = 1. / (h * h) - 1. / 3.;
            b[j] =  ((i - h) + (10. * i) + (i + h)) / 12.;
            j++;
            i += h;
        } while (j < dl - 1);

        l[j - 1] = -(fi / h);
        d[j] = (fi / h + psi);
        b[j] = -theta;
        thomasMatrix(u, d, l, dl);

        double* ptr = thomasB(u, d, l, b, dl);
        //for (int i = 0; i < dl; i++)
            //fptr2 << i * h << " " << ptr[i] << endl;
        double maxblad = fabs(ptr[0] - funkcja(0));
        for (int i = 0; i < dl; i++)
        {
            if (fabs(ptr[i] - funkcja(i * h)) > maxblad)
                maxblad = fabs(ptr[i] - funkcja(i * h));
        }
        fptr1 << log10(h) << " " << log10(maxblad) << endl;
        delete[] d;
        delete[] u;
        delete[] l;
        delete[] b;
        delete[] ptr;
        h /= 2.;
    }
    fptr1.close();
}

int main()
{
	//konwencjonalna();
    Numerow();
}
