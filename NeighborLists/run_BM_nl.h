// author luigisbailo


void run_BM_nl ( int N_A, int N_B, double R_A, double R_B, double D_A, double D_B,
				 double tau_bm, double Tsim, double L, double BMgrid ) {


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

	double boxSize = BMgrid;
	int nBoxes = (int)(L/boxSize);
	int NlimBox = (int)(10*N/nBoxes/nBoxes/nBoxes);
	if (NlimBox < 10 ) NlimBox = 10 ;
	if (NlimBox > N ) NlimBox = N ;

	struct boxcell grid [nBoxes*nBoxes*nBoxes];

	struct particle_nl particles [N];

    initPos_BM_nl ( particles, r, N_A, N_B, R_A, R_B, D_A, D_B, tau_bm, L, boxSize);

	initGridBM_nl ( particles, grid, nBoxes,  boxSize,  N);

	do {

		BFstep_nl ( particles, grid, r, nBoxes, boxSize, tau_bm, N, sqrt2TAU_BM, L );

		BFupdate_nl ( particles, N);
		
		updateGridBM_nl ( particles, grid, nBoxes, boxSize, N);

	} while ( particles[0].time < Tsim );

    gsl_rng_free (r);

}