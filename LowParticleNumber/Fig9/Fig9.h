// author luigisbailo


void fig9 (double D_A, double D_B, double R_A, double R_B) {


	double tau_bm = 0.1;
	const int N = 10; 
	const int N_A = 5;
	const int N_B = 5;

	double alpha= 9;
	double L = 25;
	int Nsamples = 10;

	int nT = 5;
	double Tsim [nT];
    Tsim[0] = 1000;
    Tsim[1] = 2000;
    Tsim[2] = 3000;
    Tsim[3] = 4000;
    Tsim[4] = 5000;

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

			run_aGF1 ( N_A, N_B, R_A, R_B, D_A, D_B, tau_bm, alpha, Tsim[t], L, stat, diffStat );

			for ( int n=0; n<N; n++){
				Diff_aGF [n][count][t] = diffStat[n];
			}

			for (int d=0; d<3; d++)
				stat[d] = 0;
			for ( int n=0; n<N; n++ )
				diffStat[n] = 0;

			run_GF ( N_A, N_B, R_A, R_B, D_A, D_B, 1., 1., tau_bm, alpha, Tsim[t], L, stat, diffStat );

			for (int n=0; n<N; n++){
				Diff_GF1 [n][count][t] = diffStat[n];
			}

			for (int d=0; d<3; d++)
				stat[d] = 0;
			for ( int n=0; n<N; n++ )
				diffStat[n] = 0;

			run_GF ( N_A, N_B, R_A, R_B, D_A, D_B, 1.5, 2.5, tau_bm, alpha, Tsim[t], L, stat, diffStat );

			for (int n=0; n<N; n++){
				Diff_GF2 [n][count][t] = diffStat[n];
			}

			for (int d=0; d<3; d++)
				stat[d] = 0;
			for ( int n=0; n<N; n++ )
				diffStat[n] = 0;

			run_BM ( N_A, N_B, R_A, R_B, D_A, D_B, tau_bm, Tsim[t], L, diffStat );

			for ( int n=0; n<N; n++){
				Diff_BM [n][count][t] = diffStat[n];
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

			}
		}
	
	sdDiff_aGF [t] = sqrt(sdDiff_aGF [t] / pow(Nsamples*N,2) );
	sdDiff_GF1 [t] = sqrt(sdDiff_GF1 [t] / pow(Nsamples*N,2) );
	sdDiff_GF2 [t] = sqrt(sdDiff_GF2 [t] / pow(Nsamples*N,2) );
	sdDiff_BM [t] = sqrt(sdDiff_BM [t] / pow(Nsamples*N,2) );

	}

	for ( int t=0; t<nT; t++){

		printf ("%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",
				Tsim[t], avDiff_aGF[t], avDiff_GF1[t], avDiff_GF2[t],
				avDiff_BM[t], sdDiff_aGF[t], sdDiff_GF1[t], sdDiff_GF2[t], sdDiff_BM[t]);

	}

}
