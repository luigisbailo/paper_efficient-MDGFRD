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
#include "../../checks.hpp"


int main (int argc, char *argv[]) {         


	high_resolution_clock::time_point tGF1,tGF2,tBM1,tBM2;
	double tGF,tBM;


	const gsl_rng_type *Type;
	gsl_rng *r;
	gsl_rng_env_setup ();
	Type = gsl_rng_default;
	r = gsl_rng_alloc (Type);
	FILE *devurandom = fopen("/dev/urandom","r");
	unsigned long int seed;
	fread(&seed, sizeof(seed), 1, devurandom);
	fclose(devurandom);
	gsl_rng_set(r, seed);

	int Nsamples;
	double dt;
	double dr; 


    stringstream convert_dt (argv[1]); 
    stringstream convert_Nsamples (argv[2]); 
    stringstream convert_dr (argv[3]); 

	if (!(convert_dt >> dt ))
		exit (EXIT_FAILURE);  
	if (!(convert_Nsamples >> Nsamples ))
		exit (EXIT_FAILURE);  
	if (!(convert_dr >> dr ))
		exit (EXIT_FAILURE);  


	double b = 0;


    double rho[10];
    for (int count=0; count<10; count++)
    	rho[count]=0;

    int idx;

    for ( int var=0; var<10; var++){

    	idx=0;

		for ( double D=0.001; D<0.011; D+=0.001) {

			b=0;
			do {

				b += dr;


		    	high_resolution_clock::time_point tGF1 = high_resolution_clock::now();

				for ( int count=0; count<Nsamples; count++) {

					 drawTimeNewt (b,D,gsl_rng_uniform(r));

				}

				high_resolution_clock::time_point tGF2 = high_resolution_clock::now();

				tGF = duration_cast<std::chrono::microseconds>(tGF2-tGF1).count();

				// cout << tGF << endl;

		    	high_resolution_clock::time_point tBM1 = high_resolution_clock::now();

				for ( int count=0; count<Nsamples; count++) {

					double dist;
					double pos[3] = {1,1,1};
					double x = 0;
					double y = 0;
					double z = 0;
					double BFvar = sqrt (2*D*dt);

					while ( x*x + y*y + z*z < b*b ){

						x += BFvar * gsl_ran_gaussian (r,1);
						y += BFvar * gsl_ran_gaussian (r,1);
						z += BFvar * gsl_ran_gaussian (r,1);

						for (int count=0; count<9; count++)
							dist = (pos[0]-x)*(pos[0]-x) + (pos[1]-y)*(pos[1]-y) + (pos[2]-z)*(pos[2]-z);



					}	


				}

				high_resolution_clock::time_point tBM2 = high_resolution_clock::now();

				tBM = duration_cast<std::chrono::microseconds>(tBM2-tBM1).count(); 
				// cout << tBM << endl;



			} while ( tGF>tBM);

			rho[idx] += b;
			idx ++;

		}
	
	}

	for (int idx=0; idx <10; idx++)
		cout << 0.001*(idx+1) << "\t" << rho[idx]/10 << endl;

}
