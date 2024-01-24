#include <cstdlib>
#include <cmath>
#include "assert.h"
#include "headerh.h"
#include <iostream>
using namespace std;


double* genz(int sz_v, double theta, unsigned int* ptr)
{
	assert(sz_v > 0);
	double rn;

	double *v = new double[sz_v];

	for(int i = 0; i < sz_v; i++)
	{
	        v[i] = 1;
		
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
