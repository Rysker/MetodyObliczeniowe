#include "Funkcje.h"
double fun1(double x, double y, double z)
{
	return x * x + y * y + z * z - 2;
}

double fun2(double x, double y, double z)
{
	return x * x + y * y - 1;
}

double fun3(double x, double y, double z)
{
	return x * x - y;
}

double df1dx1(double x)
{
	return 2 * x;
}

double df1dx2(double y)
{
	return 2 * y;
}

double df1dx3(double z)
{
	return 2 * z;
}

double df2dx1(double x)
{
	return 2 * x;
}

double df2dx2(double y)
{
	return 2 * y;
}

double df2dx3(double z)
{
	return 0;
}

double df3dx1(double x)
{
	return 2 * x;
}

double df3dx2(double y)
{
	return 0;
}

double df3dx3(double z)
{
	return 0;
}