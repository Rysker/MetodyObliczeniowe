#include <iostream>
#include <fstream>
#include <limits>
#include <string>
#include <cmath>
#include <iomanip>

using namespace std;

double maclaurin(double x)
{
    double n = 1;
    double dzielnik = 2;
    double ostatni_wyraz = 1;
    double suma = 0;
    double znak = 1;
    while(fabs(ostatni_wyraz) >= 1e-30)
    {
        suma += znak * ostatni_wyraz;
        ostatni_wyraz *=  x / dzielnik ;
        znak *= -1;
        dzielnik++;
    }
    return suma;
}


double fun(double x)
{
   return (1.0 - exp(-x))/x;
}

int main()
{
    cout.precision(53);
    cin.precision(53);
    double tab[3900][4]; // 0- log10(x) 1- x 2- (1 - exp(-x))/x 3- obliczone wyniki
    string table[3900][3];
    double wynik;
    ifstream fileptr;
    ofstream fptr;
    fileptr.open("dane.txt");
    fptr.open("wynik.txt");

    for(int i = 0; i < 3900; i++)
    {
        fileptr >> tab[i][0] >> tab[i][1] >> tab[i][2];
        double wynik = fun(tab[i][1]);

        if(tab[i][1] < 1)
            tab[i][3] = maclaurin(tab[i][1]);
        else
            tab[i][3] = wynik;
    }

    fileptr.close();
    fileptr.open("dane.txt");

     for(int i = 0; i < 3900; i++)
        fileptr >> table[i][0] >> table[i][1] >> table[i][2];

    cout << "----- Bledy -----" << endl;

    for(int i = 0; i < 3900; i++)
    {
        double error = fabs(tab[i][2] - tab[i][3]);
        if(error >= std::numeric_limits<double>::epsilon())
            cout << tab[i][3] << " " << tab[i][2] << " " << error << endl;
        fptr << setprecision(53) << table[i][0] << " " << tab[i][2] << " " << tab[i][3] << endl;
    }

    fptr.close();
    return 0;
}
