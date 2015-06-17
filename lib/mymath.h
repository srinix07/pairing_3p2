#include <iostream>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_interp.h>
#include <fstream>
#include <string>
#include <math.h>
#include <sstream>
#include <vector>
#include <gsl/gsl_linalg.h>
#include <vector>
#include <omp.h>


double interpol(double *x,double *y,int N,double x_des);
int readv(char *fi, double* &K, double** &V, int &ksize);
double rk45(double (*f)(double,double,double*),double* params, double x, double y,double max,double h);
int gauss_leg(double* x, double* w,int N,double a,double b);
double int_gauss(double (*f)(double,double*), double* params,double *x,double *w, double a,double b,int N);
double interpol2(double* x, double *V, int N,double xd,double yd);
int solv_lin(double **A, double *B, double* X, int n);
int readV(std::string fi, std::vector<double> &D);
int getV(std::vector<double> &D,std::vector<double> &K,std::vector<double> &V, int col);
int grid_split(double kf,double kstart,double kend,double offset,double offset2, int N1,int N2,int N3, double* &K, double* &W,int &N);
std::string shell_execute(std::string cmd);
