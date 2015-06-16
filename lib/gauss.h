#include <iostream>
#include <mymath.h>
#include <math.h>
#include <stdio.h>
#include <armadillo>
#include <vector>
#include <chrono>
#include <omp.h>


class gauss
{
	public:
	int size;
	std::vector<double> x,w,xr,wr;
	double a=-1,b=1;
	

	gauss();
	gauss(std::string gsize);
	void rescale(double a, double b);
	void app(gauss g);
	void print();
	std::vector<double> get(int j);
	void check();
	void stitch(gauss g1,gauss g2);
	~gauss();

};
