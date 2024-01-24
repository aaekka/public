#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include "headerh.h"
#include "omp.h"
using namespace std;

// Global Parameters
double al = 0.05;
double th = 2;
double sigma = 0.02/sqrt(2*al);
int tmax = 1000;
int h = 2;
double h2 = 1/double(h);
double pos_step = exp(sigma*sqrt(h2));
double neg_step = exp(-sigma*sqrt(h2));

int irmax = 10000;

int main(){

int N = -99;

// initialize population grid
int sz_gridN = 2;
int gridN[sz_gridN];
int N_factor = 10;
gridN[0] = 1000;
for(int i = 1; i < sz_gridN; i = i + 1)
{
	gridN[i] = gridN[i-1]*N_factor;
}

//irmax blocks of size tmax
double gr[irmax][tmax];
double gr_prv[irmax][tmax];
double gtemp[irmax];
double g[sz_gridN][tmax];

double* z;

double* intvec0;
bool* intvec;

int Ni = -99;
int Ne = -99;
int izI = -99;
double rnum;
double mean0 = -99;
double mean1 = -99;
double* Ii;
double* zI;
double start;
double end;

unsigned int uir;
unsigned int* irptr;

omp_set_num_threads(12);

// run simulations and store growth rates
for(int in = 0; in < sz_gridN; in++)
{
	N = gridN[in];

	//initialize growth rate vectors
	for(int igr = 0; igr < irmax; igr++)
	{
		for(int igr2 = 0; igr2 < tmax; igr2++)
		{
			gr[igr][igr2] = 0;
			gr_prv[igr][igr2] = 0;
		}
	}
	// avoid false-sharing by creating local copies of growth rates and aggregating at end
	start = omp_get_wtime();
	#pragma omp parallel default(none) private(z, intvec0, intvec, Ni, izI, mean0, mean1, Ii, zI, uir, irptr, rnum) firstprivate(gr_prv, N, th, irmax, tmax, al, sigma, h, pos_step, neg_step) shared(gr)
	{
		#pragma  omp for
		for(int ir = 1; ir <= irmax; ir++)
		{
			//cout << N << "  " << th << endl;
			//cout << "NumThreads = " << omp_get_num_threads() << endl;
			//cout << ir << endl;
			uir = (unsigned int) ir;
			irptr = &uir;

			intvec = new bool[N];
			
			z = genz(N,th, irptr);
	
			mean0 = vecmean(z,N);

			for(int t = 1; t <= tmax; t++)
			{
				//cout << "N = " << N << " ; ir = " << ir << " ; t = " << t << endl;

				//draw from internal distribution
				Ni = 0;
				
				intvec0 = genrandvec(N, irptr);
				
				for(int ial = 0; ial < N; ial++)
				{
					intvec[ial] = *(intvec0 + ial) < al;
					if(intvec[ial])
					{
						Ni++;
					}
				}
				
				if(Ni > 0)
				{	
					Ii = genIi(Ni, N, irptr);
					zI = new double[Ni];
					
					for(int iIi = 0; iIi < Ni; iIi++)
					{
						zI[iIi] = z[int(Ii[iIi])];
					}

					izI = 0;
					for(int iz = 0; iz < N; iz++)
					{
						if(intvec[iz])
						{
							z[iz] = fmax(z[iz], zI[izI]);
							izI++;
						}

					}
					
					delete[] Ii;
					delete[] zI;
					
				}

				delete[] intvec0;

				
				// Brownian
				for(int ih = 0; ih < h; ih++)
				{			
					for(int iz = 0; iz < N; iz++)
					{
						
						rnum = rand_r(irptr);
						rnum = fmin(rnum, RAND_MAX - 1);
						rnum = (rnum / RAND_MAX);
						if(rnum < 0.5)
						{ 
							//z[iz] = z[iz]*exp(sigma*pow(h2,0.5));
							z[iz] = z[iz]*pos_step;
						}
						else
						{ 
							//z[iz] = z[iz]*exp(-sigma*pow(h2,0.5));
							z[iz] = z[iz]*neg_step;
						}
					}			
				}
				
				
				
				// calculate growth
				mean1 = vecmean(z,N);

				gr_prv[ir-1][t-1] = (mean1 - mean0) / mean0;

				// for next time iteration
				mean0 = mean1;		
								
			}

			delete[] intvec;
			delete[] z;
		}
		
		// aggregate threads
		#pragma omp critical
		{
			for(int igr = 0; igr < irmax; igr++)
			{
				for(int igr2 = 0; igr2 < tmax; igr2++)
				{
					gr[igr][igr2] += gr_prv[igr][igr2];
				}
			}
		}	
		
	}

	end = omp_get_wtime();

	for(int it = 0; it < tmax; it++)
	{
		for(int iir = 0; iir < irmax; iir++)
		{
			gtemp[iir] = gr[iir][it];
		}
		g[in][it] = vecmean(gtemp,irmax);
		cout << "g["<<N<<"]["<<it+1<<"] = " << g[in][it] << endl;
	}

	cout << N << " done!   Time = " << (end - start) << endl;


}


//export results
ofstream outfile;
outfile.open("/home/aaekka/rslt_brn_0809small3.csv");
for(int iin = 0; iin < sz_gridN; iin++)
{
	for(int iit = 0; iit < tmax; iit++)
	{
		outfile << g[iin][iit] << ",";
	}
	outfile << endl;
}
outfile.close();


//end
}
