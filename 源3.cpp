#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>

using namespace std;

static double sigma = 10.0, beta = static_cast<double>(8.0) / 3.0, rol = 28.0;
static const double h = 0.01;
static const unsigned ITERATE_TIMES = 200000;

struct dvec3
{
	double x;
	double y;
	double z;
	dvec3(){}
	dvec3(double x, double y, double z) :
		x(x), y(y), z(z){}
};

ostream& operator<<(ostream& os, const dvec3 rhs)
{
	os << ' ' << rhs.x << ',' << rhs.y << ',' << rhs.z << ' ';
	return os;
}

dvec3 operator+(const dvec3& lhs, const dvec3& rhs)
{
	return dvec3(lhs.x + rhs.x, lhs.y + rhs.y, lhs.z + rhs.z);
}

dvec3 operator*(double factor, const dvec3& rhs)
{
	return dvec3(factor*rhs.x, factor*rhs.y, factor*rhs.z);
}

dvec3 fun(dvec3& yn)
{
	return dvec3(sigma*(yn.y - yn.x), rol*yn.x - yn.y - yn.x*yn.z, yn.x*yn.y - beta*yn.z);
}

void iterate(dvec3& yn)
{
	dvec3 k1 = fun(yn);
	dvec3 k2 = fun(yn + 0.5*h*k1);
	dvec3 k3 = fun(yn + 0.5*h*k2);
	dvec3 k4 = fun(yn + h*k3);
	yn = yn + h / static_cast<double>(6.0)*(k1 + k2 + k3 + k4);
}

int main(void)
{
	dvec3 y0 = dvec3(0, 1, 0);
	double currentT = 0.0;
	fstream fout("D:\\data.txt", ios::out);
	fout << setprecision(12);
	for (int i = 0; i != ITERATE_TIMES; ++i)
	{
		iterate(y0);
		//cout << y0 << '\t' << currentT << endl;
		fout << currentT << '\t' << y0 << endl;
		currentT += h;
	}
}