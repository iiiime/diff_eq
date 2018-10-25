// Integrate one variable ODE using Euler's method

#include <vector>
#include <string>
#include <math.h>
#include <iostream>
#include <fstream>

#define ak 0.004
#define bk 0.02
#define bs 0.082
#define k0 0.2
#define k1 0.222
#define n 2
#define p 5

#define eta 1e-7 //accuracy of error

using namespace std;

/**
 * define derivative functoin
 */
double f(double x, double y) {
	return 0;
}


double dK(double k, double s) {
	return ak + bk * pow(k, n) / (pow(k0, 2) + pow(k, n)) - k / (1 + k + s);
}

double dS(double k, double s) {
	return bs / (1 + pow((k / k1), p)) - s / (1 + k + s);
}


/**
 * x0: estimated solution
 * k: 2nd param of function dK, dS
 * 
 */
double solve(double(*f)(double, double), double x0, double x1, double k) {
	double d, cur; // difference
	while(1) {
		
		if (fabs(d) > eta && fabs(/*TODO*/) > eta) {
			break;
		}
	}
	return 0;
}

/**
 * x0, y0: initial value
 * h: step size
 */
vector<double> Euler(double x0, double y0, double h, int end){
	vector<double> y(end, 0); // initialization of the output vector
	y[0] = y0;
	double x = x0;
	for (int i = 1; i < end; i++) {
		y[i] = y[i - 1] + h * f(x, y[i - 1]);
		x += h; // increase by step size
	}
	return y;
}


int main(int argc, char* argv[]) {
	return 0;
}
