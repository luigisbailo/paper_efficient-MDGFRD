void run_aGF ( int N_A, int N_B, double R_A, double R_B, double D_A, double D_B, double tau_bm, double alpha, double Tsim, double L, double maxShell, int *stat) {

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

	particle particles [N]; 
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

	vector < vector <int> > gridGF (nBoxesGF*nBoxesGF*nBoxesGF, vector<int>(NlimBox));

 
    initPos_aGF ( particles, r, N_A, N_B, R_A, R_B, D_A, D_B, tau_bm, alpha, L, boxSizeGF ); 

    initGridGF (particles, &gridGF, NlimBox, nBoxesGF, boxSizeGF, N );

    initShell_aGF ( particles, &gridGF, distRow, distLabel, r, N, tau_bm, sqrt2TAU_BM, maxShell, L, nBoxesGF, &stat[1]);

    //sort() is a prebuild c++ funct. It sorts particles for increasing exit times
    sort ( particles, particles+N, compareTime );   
    for (int n=0; n<N; n++) partList[n]=n;

	relabel_part (particles,N);

    initGridGF (particles, &gridGF, NlimBox, nBoxesGF, boxSizeGF, N );

	int mycount = 0;
    while ( particles[partList[0]].tau_exit < Tsim ) {

		// mycount++;
		// if (mycount==10) exit(EXIT_FAILURE);

    	if ( particles[partList[0]].burst == true ) stat[0]++;


		// cout << "--------------------------------------------------------------------------------------------------------------------------\n";
		// cout << setprecision (5);
		// printPos_per ( particles, partList, N );
		// printDist_per (particles, partList, N, L);
		// cout << "\n";

    	// if ( particles[partList[0]].burst == true ) exit(EXIT_FAILURE);


		updatePart_aGF ( &particles[partList[0]], &gridGF, r, tau_bm, L, boxSizeGF, nBoxesGF );    

		// print_neighb ( particles, &gridGF,  N, nBoxesGF );


        // check_Neighb ( particles, &gridGF, N, nBoxesGF, boxSizeGF, NlimBox );


		// printPos_per ( particles, partList, N );
		// printDist_per (particles, partList, N, L);
		// cout << "\n";

		// check_aGF ( particles, partList,  N, L );

		// check_times ( particles, partList, N);

        getDist ( particles, &gridGF, distRow, distLabel, particles[partList[0]].label, N, L, nBoxesGF );

		burst_P_aGF ( particles, partList, distRow, distLabel, r, N, L); 

		R = getR_aGF ( particles, partList, distRow, distLabel, maxShell, N, L );

		particles[partList[0]].burst = false;

		if ( R > particles[partList[0]].R_bd ) {

			stat[1]++;
			if (R>L/2) R=L/2;
			GFstep_aGF ( &particles[partList[0]], r, R, L );
			particles[partList[0]].gf = true;

		}
		else{

			stat[2]++; 
			BMstep ( particles, partList, distRow, distLabel, r, tau_bm, sqrt2TAU_BM, N, L );
			particles[partList[0]].gf = false;
		}


		sortPart (particles,partList,N);

    } ;

    synchPart_P_aGF ( particles, partList, r, N, Tsim, L );

    // for ( int n=0; n<N; n++ ){

    // 	diffStat[n] += pow(particles[n].pos[0]-particles[n].pos_init[0] + particles[n].pos_period[0]*L, 2);
    // 	diffStat[n] += pow(particles[n].pos[1]-particles[n].pos_init[1] + particles[n].pos_period[1]*L, 2);
    // 	diffStat[n] += pow(particles[n].pos[2]-particles[n].pos_init[2] + particles[n].pos_period[2]*L, 2);

    // }


    gsl_rng_free (r);

}