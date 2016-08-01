#include<fstream>
#include<iostream>
#include<math.h>
#include<iomanip>
const static double PI = 3.1415926535897;
static double epsilon = 0.00001;
inline double get_fun_val(double x)
{
	return 4 / (x*x + 1);
}

double cal_rect(double h, double range_l, double range_h)
{
	int steps = (range_h - range_l) / h;
	double res = 0.0;
	for (unsigned register i = 0; i != steps; ++i)
	{
		res += get_fun_val(range_l + i*h)*h;
	}
	return res;
}

double cal_rect2(double h, double range_l, double range_h)
{
	int steps = (range_h - range_l) / h;
	double res = 0.0;
	for (unsigned register i = 1; i <= steps; ++i)
	{
		res += get_fun_val(range_l + i*h)*h;
	}
	return res;
}

double cal_stair(double h, double range_l, double range_h)
{
	int steps = (range_h - range_l) / h;
	double res = 0.0;
	for (unsigned register i = 0; i != steps; ++i)
	{
		res += (get_fun_val(range_l + i*h) + get_fun_val(range_l + (i + 1)*h))*h / 2;
	}
	return res;
}

double cal_simpson(double h, double range_l, double range_h)
{
	int steps = (range_h - range_l) / h;
	double res = 0.0;
	for (unsigned register i = 0; i != steps; ++i)
	{
		res += (get_fun_val(range_l + i*h) +4*get_fun_val(range_l+(i+(double)1/2)*h)+ get_fun_val(range_l + (i + 1)*h))*h / 6;
	}
	return res;
}

double simpson(double range_l, double range_h)
{
	double c = range_l + (range_h - range_l) / 2;
	return (get_fun_val(range_l) + 4 * get_fun_val(c) + get_fun_val(range_h))*(range_h - range_l) / 6;
}

double asr(double range_l, double range_h, double epsilon, double A)
{
	double c = range_l + (range_h - range_l) / 2;
	double L = simpson(range_l, c), R = simpson(c, range_h);
	if (fabs(A - L - R) <= 15 * epsilon) 
		return L + R + (A - L - R) / 15;
	return asr(range_l, c, epsilon / 2, L) + asr(c, range_h, epsilon / 2, R);
}
 
double asr(double range_l, double range_h, double epsilon)
{
	return asr(range_l, range_h, epsilon, get_fun_val(range_h));
}

double Romberg(double range_l, double range_h)
{
	int m, n;
	double h, x;
	double s, q;
	double ep;
	double *y = new double[10];
	double p;
	h = range_h - range_l;
	y[0] = h*(get_fun_val(range_l) + get_fun_val(range_h)) / 2.0;
	m = 1;
	n = 1;
	ep = epsilon + 1.0;
	while ((ep >= epsilon) && (m < 10))
	{
		p = 0.0;
		for (unsigned register int i = 0; i<n; i++)
		{
			x = range_l + (i + 0.5)*h;
			p = p + get_fun_val(x);
		}
		p = (y[0] + h*p) / 2.0;  
		s = 1.0;
		for (int k = 1; k <= m; k++)
		{
			s = 4.0*s;
			q = (s*p - y[k - 1]) / (s - 1.0);
			y[k - 1] = p;
			p = q;
		}
		p = fabs(q - y[m - 1]);
		m = m + 1;
		y[m - 1] = q;
		n = n + n; h = h / 2.0;
	}
	return (q);
}

int main(void)
{
	std::fstream fout("D:\\data_rect.txt", std::ios::out);
	fout << std::setprecision(14);
	
	epsilon = 1e-9;
	while (epsilon < 1e-1)
	{
		fout << epsilon << ',' << asr(0.0, 1.0,epsilon) - PI << '\n';
		epsilon *= 5;
		fout << epsilon << ',' << asr(0.0, 1.0, epsilon) - PI << '\n';
		epsilon *= 2;
		fout << epsilon << ',' << asr(0.0, 1.0, epsilon) - PI << '\n';
	}
	fout.close();
	std::cout << asr(0.0,1.0,0.01);
}