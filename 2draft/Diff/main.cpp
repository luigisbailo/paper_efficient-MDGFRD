//TO COMPILE: g++ -std=c++11 main.cpp -o main -lgsl -lgslcblas -lm

using namespace std;
#define print(x) cout << x << endl;

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_math.h>
#include <cmath>
#include <gsl/gsl_errno.h>
#include <iostream>
#include <fstream>
#include <algorithm> 
#include <vector>
#include <iomanip>
#include <chrono>
#include <cstring>
#include <sstream>

using namespace std;
using namespace std::chrono;


#include "../parameters.hpp"
#include "../tools.hpp"
#include "../greensFunct.hpp"
#include "../draw.hpp"
#include "../init.hpp"
#include "../step.hpp"
#include "../shell.hpp"
#include "../print.hpp"
#include "../burst.hpp"
#include "../bruteForce.hpp"
#include "../checks.hpp"
#include "../run_aGF.hpp"
#include "../run_GF.hpp"
#include "../run_BM.hpp"


int main (int argc, char *argv[]) {


	
	double D_A;
	double D_B;
	double R_A;
	double R_B;
	double MU_BM;
	double MU_GF;
	
	double tau_bm = 0.1;
	const int N = 10; 
	const int N_A = 5;
	const int N_B = 5;

	double alpha= 9;
	double L = 25;
	// int Nsamples = 100;
	int Nsamples = 1000;

	int nT = 5;
	double Tsim [nT] = {100000,200000,300000,400000,500000};
	// int nT = 1;
	// double Tsim [nT] = {1000};

	stringstream convert_DA (argv[1]); 
	stringstream convert_RA (argv[2]); 
	stringstream convert_DB (argv[3]); 
	stringstream convert_RB (argv[4]); 

	if (!(convert_DA >> D_A ))
	exit (EXIT_FAILURE);  
	if (!(convert_RA >> R_A ))
	exit (EXIT_FAILURE);  
	if (!(convert_DB >> D_B ))
	exit (EXIT_FAILURE);  
	if (!(convert_RB >> R_B ))
	exit (EXIT_FAILURE);  

	int stat[3];
	double diffStat [N];
	double Diff_aGF [N][Nsamples][nT];
	double Diff_GF1  [N][Nsamples][nT];
	double Diff_GF2  [N][Nsamples][nT];
	int Diff_BM  [N][Nsamples][nT];

	for ( int t=0; t<nT; t++) {

		for ( int count = 0; count < Nsamples; count++){


			for (int d=0; d<3; d++ )
				stat[d] = 0;
			for ( int n=0; n<N; n++ )
				diffStat[n] = 0;

			run_aGF ( N_A, N_B, R_A, R_B, D_A, D_B, tau_bm, alpha, Tsim[t], L, stat, diffStat );

			for ( int n=0; n<N; n++){
				Diff_aGF [n][count][t] = diffStat[n];
				// cout << Diff_aGF [n][count][t]<<"\t";
			}



			for (int d=0; d<3; d++)
				stat[d] = 0;
			for ( int n=0; n<N; n++ )
				diffStat[n] = 0;

			run_GF ( N_A, N_B, R_A, R_B, D_A, D_B, 1., 1., tau_bm, alpha, Tsim[t], L, stat, diffStat );

			for (int n=0; n<N; n++){
				Diff_GF1 [n][count][t] = diffStat[n];
				// cout << Diff_GF [n][count][t]<<"\t";
			}



			for (int d=0; d<3; d++)
				stat[d] = 0;
			for ( int n=0; n<N; n++ )
				diffStat[n] = 0;

			run_GF ( N_A, N_B, R_A, R_B, D_A, D_B, 1.5, 2.5, tau_bm, alpha, Tsim[t], L, stat, diffStat );

			for (int n=0; n<N; n++){
				Diff_GF2 [n][count][t] = diffStat[n];
				// cout << Diff_GF [n][count][t]<<"\t";
			}



			for (int d=0; d<3; d++)
				stat[d] = 0;
			for ( int n=0; n<N; n++ )
				diffStat[n] = 0;

			run_BM ( N_A, N_B, R_A, R_B, D_A, D_B, tau_bm, Tsim[t], L, diffStat );

			for ( int n=0; n<N; n++){
				Diff_BM [n][count][t] = diffStat[n];
				// cout << Diff_BM [n][count][t]<<"\t";
			}
		}	
	
	}	

	double avDiff_aGF [nT];
	double avDiff_GF1 [nT];
	double avDiff_GF2 [nT];
	double avDiff_BM [nT];

	for ( int t=0; t<nT; t++ ){

		avDiff_aGF [t] = 0;
		avDiff_GF1 [t] = 0;
		avDiff_GF2 [t] = 0;
		avDiff_BM [t] = 0;

		for ( int count=0; count<Nsamples; count++){
			for ( int n=0; n<N; n++) {

				avDiff_aGF [t] += Diff_aGF [n][count][t]; 
				avDiff_GF1 [t] += Diff_GF1 [n][count][t]; 
				avDiff_GF2 [t] += Diff_GF2 [n][count][t]; 
				avDiff_BM [t] += Diff_BM [n][count][t]; 
			}
		}
	
	avDiff_aGF [t] = avDiff_aGF [t] / (Nsamples*N);
	avDiff_GF1 [t] = avDiff_GF1 [t] / (Nsamples*N);
	avDiff_GF2 [t] = avDiff_GF2 [t] / (Nsamples*N);
	avDiff_BM [t] = avDiff_BM [t] / (Nsamples*N);

	}


	double sdDiff_aGF [nT];
	double sdDiff_GF1 [nT];
	double sdDiff_GF2 [nT];
	double sdDiff_BM [nT];

	for ( int t=0; t<nT; t++ ){

		sdDiff_aGF [t] = 0;
		sdDiff_GF1 [t] = 0;
		sdDiff_GF2 [t] = 0;
		sdDiff_BM [t] = 0;

		for ( int count=0; count<Nsamples; count++){
			for (int n=0; n<N; n++){

				sdDiff_aGF [t] += pow(Diff_aGF[n][count][t]-avDiff_aGF[t],2); 
				sdDiff_GF1 [t] += pow(Diff_GF1[n][count][t]-avDiff_GF1[t],2); 
				sdDiff_GF2 [t] += pow(Diff_GF2[n][count][t]-avDiff_GF2[t],2); 
				sdDiff_BM [t] += pow(Diff_BM[n][count][t]-avDiff_BM[t],2); 

			
				// cout <<  pow(Diff_aGF[n][count][t]-avDiff_aGF[t],2) 
				// 	<<"\t"<< pow(Diff_GF[n][count][t]-avDiff_GF[t],2) 
				// 	<< "\t"<<pow(Diff_BM[n][count][t]-avDiff_BM[t],2)<<endl; 
			
			}
		}
	
	sdDiff_aGF [t] = sqrt(sdDiff_aGF [t] / pow(Nsamples*N,2) );
	sdDiff_GF1 [t] = sqrt(sdDiff_GF1 [t] / pow(Nsamples*N,2) );
	sdDiff_GF2 [t] = sqrt(sdDiff_GF2 [t] / pow(Nsamples*N,2) );
	sdDiff_BM [t] = sqrt(sdDiff_BM [t] / pow(Nsamples*N,2) );

	}



	cout << setprecision (7);

	for ( int t=0; t<nT; t++){

		cout << Tsim[t] << "\t" << avDiff_aGF[t] << "\t" << avDiff_GF1[t] << "\t" << avDiff_GF2[t] << "\t" << avDiff_BM[t] << "\t" ;
		cout << sdDiff_aGF[t] << "\t" << sdDiff_GF1[t] << "\t" << sdDiff_GF2[t] << "\t" << sdDiff_BM[t] << endl;

	}

}
