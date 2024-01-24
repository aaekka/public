#ifndef HEADERH_H
#define HEADERH_H

double* genz(int sz_v, double theta, unsigned int* ptr, double eps, double L);

double* genlnz(int sz_v, double theta, unsigned int* ptr, double mu, double sigma);

double* genrandvec(int sz_v, unsigned int* ptr);

double* genIi(int sz_v, int N, unsigned int* ptr);

double* genztemp(int sz_v, double theta, unsigned int* ptr);

double vecmean(double* v, int sz_v);

double vecstd(double* v, int sz_v);

#endif