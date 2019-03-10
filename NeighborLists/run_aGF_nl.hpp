void run_aGF_nl ( int N_A, int N_B, double R_A, double R_B, double D_A, double D_B, double tau_bm,
				  double alpha, double Tsim, double L, double maxShell, int *stat) {

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

	double R;
	int nBurst = 0;
	int nGF = 0;

	const double sqrt2TAU_BM = sqrt(2*tau_bm);

	const int N = N_A + N_B;

	particle_nl particles [N];
	double distRow [N];
	int distLabel [N];
	for (int n=0; n<N; n++) {
		distRow[n]=-1;
		distLabel[n]=-1;
	}
	int partList [N];

	double boxSizeGF = (R_A+maxShell)*2;

	int nBoxesGF = int(L/boxSizeGF);
	int NlimBox = int(10*N/nBoxesGF/nBoxesGF/nBoxesGF);
	if (NlimBox < 10 ) NlimBox = 10 ;
	if (NlimBox > N ) NlimBox = N ;

	std::vector < std::vector <int> > gridGF (nBoxesGF*nBoxesGF*nBoxesGF, std::vector<int>(NlimBox));

    initPos_aGF_nl ( particles, r, N_A, N_B, R_A, R_B, D_A, D_B, tau_bm, alpha, L, boxSizeGF );

    initGridGF_nl (particles, &gridGF, NlimBox, nBoxesGF, boxSizeGF, N );

    initShell_aGF_nl ( particles, &gridGF, distRow, distLabel, r, N, tau_bm, sqrt2TAU_BM, maxShell, L, nBoxesGF, &stat[1]);

    std::sort ( particles, particles+N, compareTime_nl );
    for (int n=0; n<N; n++) partList[n]=n;

	relabel_part_nl (particles,N);

    initGridGF_nl (particles, &gridGF, NlimBox, nBoxesGF, boxSizeGF, N );

	int mycount = 0;
    while ( particles[partList[0]].tau_exit < Tsim ) {

    	if ( particles[partList[0]].burst == true ) stat[0]++;

		updatePart_aGF_nl ( &particles[partList[0]], &gridGF, r, tau_bm, L, boxSizeGF, nBoxesGF );

        getDist_nl ( particles, &gridGF, distRow, distLabel, particles[partList[0]].label, N, L, nBoxesGF );

		burst_P_aGF_nl ( particles, partList, distRow, distLabel, r, N, L);

		R = getR_aGF_nl ( particles, partList, distRow, distLabel, maxShell, N, L );

		particles[partList[0]].burst = false;

		if ( R > particles[partList[0]].R_bd ) {

			stat[1]++;
			if (R>L/2) R=L/2;
			GFstep_aGF_nl ( &particles[partList[0]], r, R, L );
			particles[partList[0]].gf = true;

		}
		else{

			stat[2]++; 
			BMstep_nl ( particles, partList, distRow, distLabel, r, tau_bm, sqrt2TAU_BM, N, L );
			particles[partList[0]].gf = false;
		}


		sortPart_nl (particles,partList,N);

    } ;

    synchPart_P_aGF_nl ( particles, partList, r, N, Tsim, L );

    gsl_rng_free (r);

}
