void run_BM ( int N_A, int N_B, double R_A, double R_B, double D_A, double D_B, double tau_bm, double Tsim, double L, double BMgrid ) {


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

	// double boxSize = R_A*2;
	double boxSize = BMgrid;
	int nBoxes = int(L/boxSize);
	int NlimBox = int(10*N/nBoxes/nBoxes/nBoxes);
	if (NlimBox < 10 ) NlimBox = 10 ;
	if (NlimBox > N ) NlimBox = N ;

	vector < vector <int> > grid (nBoxes*nBoxes*nBoxes, vector<int>(NlimBox));

	// int partList [N];
	// for (int n=0; n<N; n++)partList[n]=n;

	// BFdistances d[N];
	particle particles [N]; 
	// cout << D_A << "\t" << D_B << endl;
   

    initPos_BM ( particles, r, N_A, N_B, R_A, R_B, D_A, D_B, tau_bm, L, boxSize); 

	initGridBM ( particles, &grid, NlimBox, nBoxes,  boxSize,  N);

	check_Neighb ( particles, &grid, N, nBoxes, boxSize, NlimBox );


	// int mycount = 0;

	do {
		// mycount ++;
		// if (mycount==10) exit (EXIT_FAILURE);

		BFstep ( particles, &grid, r, nBoxes, boxSize, tau_bm, N, sqrt2TAU_BM, L );


		// cout << "--------------------------------------------------------------------------------------------------------------------------\n";
		// cout << setprecision(6);
		// printPos_per ( particles, partList, N );
		// printDist_per (particles, partList, N, L);
		// cout << "\n";


		BFupdate ( particles, N);
		
		updateGridBM ( particles, &grid, nBoxes, boxSize, N);

	} while ( particles[0].time < Tsim );


    // for ( int n=0; n<N; n++ ){

    // 	diffStat[n] += pow(particles[n].pos[0]-particles[n].pos_init[0] + particles[n].pos_period[0]*L, 2);
    // 	diffStat[n] += pow(particles[n].pos[1]-particles[n].pos_init[1] + particles[n].pos_period[1]*L, 2);
    // 	diffStat[n] += pow(particles[n].pos[2]-particles[n].pos_init[2] + particles[n].pos_period[2]*L, 2);

    // }


    gsl_rng_free (r);


}





