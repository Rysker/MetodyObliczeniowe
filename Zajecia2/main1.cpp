#include <iostream>
#include <fstream>

using namespace std;
double fun(double x)
{
    return (1.0-exp(-x))/x;
}


int main()
{
    double tab[3900][4]; // 0- log10(x) 1- x 2- (1 - exp(-x))/x 3- obliczone wyniki
    double wynik;
    ifstream fileptr;
    ofstream fptr;
    fileptr.open("dane.txt");
    fptr.open("wynik.txt");
    for(int i = 0; i < 3900; i++)
    {
        fileptr >> tab[i][0] >> tab[i][1] >> tab[i][2];
        tab[i][3] = fun(tab[i][1]);
    }

    for(int i = 0; i < 3900; i++)
    {
        cout << tab[i][2] << " " << tab[i][1] << " " << tab[i][3] << endl;
        double error = abs((tab[i][3] - tab[i][2])/tab[i][2]);
        error = log10(error);
        fptr << tab[i][0] << " " << error << endl;
    }
    return 0;
}
