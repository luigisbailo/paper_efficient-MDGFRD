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


#include "../../parameters.hpp"
#include "../../tools.hpp"
#include "../../greensFunct.hpp"
#include "../../draw.hpp"
#include "../../init.hpp"
#include "../../step.hpp"
#include "../../shell.hpp"
#include "../../print.hpp"
#include "../../burst.hpp"
#include "../../bruteForce.hpp"
#include "../../run_aGF2.hpp"
#include "../../run_GF.hpp"


int main ( int argc, char *argv[] ) {


	double L;
	double D_A;
	double D_B;
	const double R_A = 2.5;
	const double R_B = 2.5;
	double tau_bm;
	int nINTsteps;
	int POWsteps;
	int Nsamples;
	const int N = 10;	
	const int N_A = 5;
	const int N_B = 5;

	int nAlphas = 6;
	double alphaValues[nAlphas] = {6,8,10,12,14,16};

	int stat [3];
	double diffStat [N];

    stringstream convert_L (argv[1]); 
    stringstream convert_DA (argv[2]); 
    stringstream convert_DB (argv[3]); 
    stringstream convert_TAUBM (argv[4]); 
    stringstream convert_POWSTEPS (argv[5]); 
    stringstream convert_SAMPLES (argv[6]); 

	if (!(convert_L >> L ))
		exit (EXIT_FAILURE);  
	if (!(convert_DA >> D_A ))
		exit (EXIT_FAILURE);  
	if (!(convert_DB >> D_B ))
		exit (EXIT_FAILURE);  
	if (!(convert_TAUBM >> tau_bm ))
		exit (EXIT_FAILURE);  
	if (!(convert_POWSTEPS >> POWsteps ))
		exit (EXIT_FAILURE);  
	if (!(convert_SAMPLES >> Nsamples ))
		exit (EXIT_FAILURE);  


	nINTsteps = pow (10,POWsteps);
	const double Tsim = nINTsteps*tau_bm;


	double arrStat [3][nAlphas];	
	double arrT [nAlphas];	
	// double arrT [nAlphas][Nsamples];	
	// double avTs[nAlphas], varTs[nAlphas];

	for ( int count=0; count<nAlphas; count++){
		arrT [count] = 0;
		// avTs [count] = 0;
		// varTs [count] = 0;
	}
	for (int d=0; d<3; d++) {
		for (int count=0; count<nAlphas; count++) {
			arrStat[d][count]=0;
		}
	}


	high_resolution_clock::time_point t2,t1;
	double t12, avT;


	cout << D_A << endl;
    cout << D_B << endl;
    cout << L << endl;
    cout << tau_bm << endl;
    cout << nINTsteps << endl;
    cout << Nsamples << endl;
    

	for ( int count=0; count < Nsamples; count ++){

		for ( int  n=0; n<nAlphas; n++ ){
 
			double alpha = alphaValues[n];

 			for (int d=0; d<3; d++) 
 				stat [d] = 0;

 			for (int n=0; n<N; n++)
 				diffStat [N] = 0;

		    t1 = high_resolution_clock::now();

		    run_aGF2 ( N_A, N_B, R_A, R_B, D_A, D_B, tau_bm, alpha, Tsim, L, stat, diffStat );

		    t2 = high_resolution_clock::now();
		    
		    t12 = double(std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count())/1000000;

		    arrT [n] += t12;
		    arrStat [0][n] += stat[0]; 
		    arrStat [1][n] += stat[1]; 
		    arrStat [2][n] += stat[2]; 

		    // stat[0] is the number of domains burst
		    // stat[1] is the number of domains constructed
		    // stat[2] is the number of brute-force stepsls

		}

	}

	cout << setprecision (9);

	// for (int n=0; n<nAlphas; n++) {
	// 	for (int count=0; count<Nsamples; count++) {
	// 		avTs [n] += arrT [n][count];
	// 	}
	// 	avTs [n] = avTs[n]/Nsamples;
	// }


	// for (int n=0; n<nAlphas; n++) {
	// 	for (int count=0; count<Nsamples; count++) {
	// 		varTs [n] += pow(arrT [n][count] - avTs[n],2);
	// 	}
	// 	varTs [n] = sqrt(vaTrs[n])/Nsamples;
	// }



	for (int count=0; count<nAlphas; count++)
		cout << alphaValues[count] << "\t" << arrT[count]/Nsamples  << "\t" 
			 << double(arrStat[0][count])/Nsamples << "\t" << double(arrStat[1][count])/Nsamples << "\t" << double(arrStat[2][count])/Nsamples << endl;

 
}
