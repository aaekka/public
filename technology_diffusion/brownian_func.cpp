#include <cstdlib>
#include <cmath>
#include "assert.h"
#include "headerh.h"
#include <iostream>
using namespace std;


double* genz(int sz_v, double theta, unsigned int* ptr, double eps, double L)
{
	assert(sz_v > 0);

	double *v = new double[sz_v];

	for(int i = 0; i < sz_v; i++)
	{
		if (i < eps*sz_v) {
        	v[i] = exp(L);
		}
		else{
			v[i] = 1;
		}
	}

	return v;
}

double* genlnz(int sz_v, double theta, unsigned int* ptr, double mu, double sigma)
{
	assert(sz_v > 0);

	double *v = new double[sz_v];
	double u1;
	double u2;
	double z1;
	double z2;
	bool newdraw = true;
	const double pi = 3.14159265358979323846;

	for(int i = 0; i < sz_v; i++)
	{
		// use Box-Muller transformation to get lognormal variables
		if (newdraw == true) {
			// get uniform RV
        	u1 = rand_r(ptr);
			u2 = rand_r(ptr);
			u1 = fmin(u1, RAND_MAX - 1);
			u2 = fmin(u2, RAND_MAX - 1);
			u1 = (u1 / RAND_MAX);
			u2 = (u2 / RAND_MAX);

			// convert to standard normal RV
			z1 = sqrt(-2*log(u1))*cos(2*pi*u2);
			z2 = sqrt(-2*log(u1))*sin(2*pi*u2);
			// convert to normal RV with proper sigma and mu
			z1 = (z1 * sigma) + mu;
			z2 = (z2 * sigma) + mu;

			v[i] = exp(z1);
			newdraw = false;
		}
		else{
			v[i] = exp(z2);
			newdraw = true;
		}
	}

	return v;
}

double* genrandvec(int sz_v, unsigned int* ptr)
{
	assert(sz_v > 0);
	double rn;

	double * v = new double[sz_v];

	for(int i = 0; i < sz_v; i++)
	{
		//to ensure random number is not 1
		rn = rand_r(ptr);
		rn = fmin(rn, RAND_MAX - 1);
		rn = (rn / RAND_MAX);

		v[i] = rn;
	}

	return v;

}

double* genIi(int sz_v, int N, unsigned int* ptr)
{
	assert(sz_v > 0);
	double rn;

	double * v = new double[sz_v];

	for(int i = 0; i < sz_v; i++)
	{
		//to ensure random number is not 1
		rn = rand_r(ptr);
		rn = fmin(rn, RAND_MAX - 1);
		rn = (rn / RAND_MAX);

		v[i] = fmax(0,fmin(N,round(N*rn)));
	}

	return v;
}

double* genztemp(int sz_v, double theta, unsigned int* ptr)
{
	assert(sz_v > 0);
	double rn;

	double * v = new double[sz_v];

	for(int i = 0; i < sz_v; i++)
	{
		//to ensure random number is not 1
		rn = rand_r(ptr);
		rn = fmin(rn, RAND_MAX - 1);
		rn = (rn / RAND_MAX);

		v[i]  = pow((1/(1-rn)),(1/theta));
	}

	return v;
}


double vecmean(double* v, int sz_v)
{
	assert(sz_v > 0);
	
	double sum = 0;
	
	for(int i = 0; i < sz_v; i++)
	{
		sum = sum + *v;
		v++;
	}

	return sum/sz_v;
}

double vecstd(double* v, int sz_v)
{
	assert(sz_v > 0);
	double sse = 0;

	double mean = vecmean(v, sz_v);

	for(int i = 0; i < sz_v; i++)
	{
		sse = sse + pow(*v - mean,2);
	}

	return sqrt(sse/sz_v);
}
