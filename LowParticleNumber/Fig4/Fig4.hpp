// author luigisbailo


void fig4 () {


	std::chrono::high_resolution_clock::time_point tGF1,tGF2,tBM1,tBM2;
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

	int Nsamples=100;
	double dt=0.001;
	double dr=0.001;
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


		    	tGF1 = std::chrono::high_resolution_clock::now();

				for ( int count=0; count<Nsamples; count++) {

					 drawTimeNewt (b,D,gsl_rng_uniform(r));

				}

				tGF2 = std::chrono::high_resolution_clock::now();

				tGF = std::chrono::duration_cast<std::chrono::microseconds>(tGF2-tGF1).count();

				// cout << tGF << endl;

		    	tBM1 = std::chrono::high_resolution_clock::now();

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

				tBM2 = std::chrono::high_resolution_clock::now();

				tBM = std::chrono::duration_cast<std::chrono::microseconds>(tBM2-tBM1).count();
				// cout << tBM << endl;



			} while ( tGF>tBM);

			rho[idx] += b;
			idx ++;

		}
	
	}

	for (int idx=0; idx <10; idx++)
		std::cout << 0.001*(idx+1) << "\t" << rho[idx]/10 << std::endl;

}
