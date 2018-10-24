// Integrate one variable ODE using Euler's method

#include <vector>
#include <string>
#include <iostream>
#include <fstream>

using namespace std;

/**
 * define derivative functoin
 */
double f(double x, double y){
	//TODO
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

