// author luigisbailo


void run_BM ( int N_A, int N_B, int R_A, int R_B, double D_A, double D_B, double tau_bm, double Tsim, double L, double *diffStat ) {


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

	const double sqrt2TAU_BM = sqrt(2*tau_bm);

	const int N = N_A + N_B;

	struct BFdistances d[N];
	struct particle particles [N];

    initPos_BM ( particles, r, N_A, N_B, R_A, R_B, D_A, D_B, tau_bm, L); 

	do {

		getBFdistances ( particles, d, N ,L );

    	BFstep ( particles, d, r, tau_bm, N, sqrt2TAU_BM, L );

		BFupdate ( particles, N);

	} while ( particles[0].time < Tsim );

    for ( int n=0; n<N; n++ ){

    	diffStat[n] += pow(particles[n].pos[0]-particles[n].pos_init[0] + particles[n].pos_period[0]*L, 2);
    	diffStat[n] += pow(particles[n].pos[1]-particles[n].pos_init[1] + particles[n].pos_period[1]*L, 2);
    	diffStat[n] += pow(particles[n].pos[2]-particles[n].pos_init[2] + particles[n].pos_period[2]*L, 2);

    }

    gsl_rng_free (r);

}





