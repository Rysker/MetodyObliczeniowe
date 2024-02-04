#include <iostream>

using namespace std;

int main()
{
    int bity = 0;
    float epsilon_float = 1.0f;
    double epsilon_double = 1.0;

    while((1.f + epsilon_float/2.f) > 1.f)
    {
        epsilon_float /= 2.f;
        bity++;
    }

    cout << "Float\n" << "Epsilon: " << epsilon_float << "\nBity: " << bity << endl << endl;

    bity = 0;

    while((1.0 + epsilon_double/2.0) > 1.0)
    {
        epsilon_double /= 2.0;
        bity++;
    }

    cout << "Double\n" << "Epsilon: " << epsilon_double << "\nBity: " << bity << endl;
    return 0;
}
