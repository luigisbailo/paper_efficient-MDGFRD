void fig4 () {


	clock_t start_tGF, end_tGF, total_tGF, start_tBM, end_tBM, total_tBM;

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


				start_tGF = clock();

				for ( int count=0; count<Nsamples; count++) {

					 drawTimeNewt (b,D,gsl_rng_uniform(r));

				}

				end_tGF = clock();

				total_tGF = (double)(end_tGF - start_tGF) ;

				start_tBM = clock();

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

						for (int count2=0; count2<9; count2++)
							dist = (pos[0]-x)*(pos[0]-x) + (pos[1]-y)*(pos[1]-y) + (pos[2]-z)*(pos[2]-z);



					}	


				}

				end_tBM = clock();

				total_tBM = (double)(end_tBM - start_tBM);



			} while ( total_tGF>total_tBM);

			rho[idx] += b;
			idx ++;

		}
	
	}

	for (int idx=0; idx <10; idx++) {

		printf("%lf\t%lf\n", 0.001 * (idx + 1), rho[idx] / 10);

	}

}
