// Integrate one variable ODE using Euler's method

#include <vector>
#include <string>
#include <math.h>
#include <iostream>
#include <fstream>

#define ak 0.004
#define bk 0.07
#define bs 0.82
#define k0 0.2
#define k1 0.222
#define n 2
#define p 5

#define eta 1e-9 //accuracy of error

using namespace std;

/**
 * define derivative functoin
 */
double f(double x, double y) {
	return 0;
}


double dK(double s, double k) {
	return ak + bk * pow(k, n) / (pow(k0, n) + pow(k, n)) - k / (1 + k + s);
}

double dS(double s, double k) {
	return bs / (1 + pow((k / k1), p)) - s / (1 + k + s);
}

double test(double k, double s) {
	return s*s + k;
}
/**
 * [x0, x1]: estimated interval of solution
 * k: param of function dK, dS
 * 
 */
double solve(double(*f)(double, double), double x1, double x0, double k) {
	double d, x2; // difference
	int i = 0;
	while(1) {
		x2 = x1 - (*f)(x1, k) * (x1 - x0) / ((*f)(x1, k) - (*f)(x0, k));
		d = x1 - x0;
		//cout << "x2: " << x2 << " " << "x1: " << x1 << endl;
		//cout << d << endl;
		if (fabs(d) < eta && fabs((*f)(x1, k)) < eta) {
			break;
		}
		x0 = x1;
		x1 = x2;
	}
	return x1;
}

/**
 * x0, y0: initial value
 * h: step size
 */
vector<double> Euler(double x0, double y0, double h, int end){
	vector<double> y(end, 0); // initialization of the output vector
	y[0] = y0;
	double x = x0;
	for (double i = 1; i < end; i++) {
		y[i] = y[i - 1] + h * f(x, y[i - 1]);
		x += h; // increase by step size
	}
	return y;
}


int main(int argc, char* argv[]) {
	for(int i = 0 ; i < 450; i ++) {
		double x = i * 0.001;
		double y_s = solve(dS, 0, 6, x);
		double y_k = solve(dK, 0, 6, x);
		cout << x << " " << y_k << " " << y_s << endl;
	}
	//double x = solve(dK, 0, 6, 0.2);
	//cout << x << endl;
	return 0;
}
