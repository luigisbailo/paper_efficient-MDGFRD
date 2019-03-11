// author luigisbailo


void fig5_8 ( double D_A, double D_B, double R_A, double R_B ) {

    int nL = 7;
    double L[nL];
    L[0]=2250;
    L[1]=1150;
    L[2]=545;
    L[3]=255;
    L[4]=118;
    L[5]=55;
    L[6]=25;
    int Nsamples[nL];
//    Nsamples[0]=1000;
//    Nsamples[1]=1000;
//    Nsamples[2]=1000;
//    Nsamples[3]=100;
//    Nsamples[4]=100;
//    Nsamples[5]=10;
//    Nsamples[6]=1;

    Nsamples[0]=100;
    Nsamples[1]=100;
    Nsamples[2]=100;
    Nsamples[3]=10;
    Nsamples[4]=10;
    Nsamples[5]=1;
    Nsamples[6]=1;


    double tau_bm=0.1;
    long int nINTsteps;
    int POWsteps=7;
    const int N = 10;
    const int N_A = 5;
    const int N_B = 5;
    double alpha = 9.;
    int BMsamples = 1;


    clock_t start_t, end_t, total_t;
    double t12, avT;

    nINTsteps = (int) (pow (10,POWsteps));
    const double Tsim = nINTsteps*tau_bm;

    int stat [3];
    double diffStat[N];

    // stat[0] number burstings
    // stat[1] number domains
    // stat[2] BD steps
    // stat[3] squared displacement
    double stat_aGF1[3][nL][Nsamples[0]*BMsamples];
    double stat_aGF2[3][nL][Nsamples[0]*BMsamples];
    double stat_hybGF[3][nL][Nsamples[0]*BMsamples];
    double stat_GF1[3][nL][Nsamples[0]*BMsamples];
    double stat_GF2[3][nL][Nsamples[0]*BMsamples];
    double stat_BM [nL][BMsamples];
    double arrT [6][nL][Nsamples[0]*BMsamples];

    printf ("%lf\n%lf\n%lf\n%ld\n%lf\n",D_A,D_B,tau_bm,nINTsteps,alpha);

    for ( int count1=0; count1 < BMsamples; count1 ++){

        int n=0;

        for ( int l=0; l<nL; l++) {

            for ( int count2=0; count2<Nsamples[l]; count2++ ) {

                for (int d=0; d<3; d++)
                    stat[d] = 0;

                start_t = clock();

                run_aGF1 ( N_A, N_B, R_A, R_B, D_A, D_B, tau_bm, alpha, Tsim, L[l], stat, diffStat );

                end_t = clock();

                total_t = (double)(end_t - start_t) / CLOCKS_PER_SEC;

                arrT [0][l][count1*Nsamples[l]+count2] = total_t;

                for (int d=0; d<3; d++)
                    stat_aGF1 [d][l][count1*Nsamples[l]+count2] = stat[d];

                for (int d=0; d<3; d++)
                    stat[d] = 0;

                for (int n=0; n<N; n++)
                    diffStat[n]= 0;

                start_t = clock();

//                run_aGF2 ( N_A, N_B, R_A, R_B, D_A, D_B, tau_bm, alpha, Tsim, L[l], stat, diffStat );

                end_t = clock();

                total_t = (double)(end_t - start_t) / CLOCKS_PER_SEC;

                arrT [1][l][count1*Nsamples[l]+count2] = total_t;

                for (int d=0; d<3; d++)
                    stat_aGF2 [d][l][count1*Nsamples[l]+count2] = stat[d];

                for (int d=0; d<3; d++)
                 stat[d] = 0;
                for (int n=0; n<N; n++)
                 diffStat[n]=0;

                start_t = clock();

//                run_hybGF ( N_A, N_B, R_A, R_B, D_A, D_B, tau_bm, alpha, Tsim, L[l], stat, diffStat );

                end_t = clock();

                total_t = (double)(end_t - start_t) / CLOCKS_PER_SEC;

                arrT [2][l][count1*Nsamples[l]+count2] = total_t;

                for (int d=0; d<3; d++)
                    stat_hybGF [d][l][count1*Nsamples[l]+count2] = stat[d];

                for (int d=0; d<3; d++)
                  stat[d] = 0;
                for (int n=0; n<N; n++)
                  diffStat[n] = 0;

                start_t = clock();

//                run_GF ( N_A, N_B, R_A, R_B, D_A, D_B, 1.,1., tau_bm, alpha, Tsim, L[l], stat, diffStat );

                end_t = clock();

                total_t = (double)(end_t - start_t) / CLOCKS_PER_SEC;


                arrT [3][l][count1*Nsamples[l]+count2] = total_t;

                for (int d=0; d<3; d++)
                  stat_GF1 [d][l][count1*Nsamples[l]+count2] = stat[d];


                for (int d=0; d<3; d++)
                    stat[d] = 0;
                for (int n=0; n<N; n++)
                diffStat[n] = 0;

                start_t = clock();

//                run_GF ( N_A, N_B, R_A, R_B, D_A, D_B, 1.5, 2.5, tau_bm, alpha, Tsim, L[l], stat, diffStat );

                end_t = clock();

                total_t = (double)(end_t - start_t) / CLOCKS_PER_SEC;

                arrT [4][l][count1*Nsamples[l]+count2] = total_t;

                for (int d=0; d<3; d++)
                    stat_GF2 [d][l][count1*Nsamples[l]+count2] = stat[d];
            }


        for (int n=0; n<N; n++)
        diffStat[n] = 0;

        start_t = clock();

//        run_BM ( N_A, N_B, R_A, R_B, D_A, D_B, tau_bm, Tsim, L[l], diffStat );

        end_t = clock();


        arrT [5][l][count1] = t12;

//        stat_BM [l][count1] = stat[3];

        }


    }

  // averages computing
  double stat_aGF1av[3][nL],stat_aGF2av[3][nL],stat_hybGFav[3][nL], stat_GF1av[3][nL],stat_GF2av[3][nL],stat_BMav[nL];
  for ( int l=0; l<nL; l++) {
    for ( int d=0; d<3; d++ ){
      stat_GF1av[d][l] = 0;
      stat_GF2av[d][l] = 0;
      stat_aGF1av [d][l] = 0;
      stat_aGF2av [d][l] = 0;
      stat_hybGFav [d][l] = 0;
    }
    stat_BMav [l] =0;
  }
  for ( int l=0; l<nL; l++) {
    for ( int d=0; d<3; d++){

      for ( int count=0; count<Nsamples[l]*BMsamples; count++){
        stat_aGF1av [d][l] += stat_aGF1 [d][l][count];
        stat_aGF2av [d][l] += stat_aGF2 [d][l][count];
        stat_hybGFav [d][l] += stat_hybGF [d][l][count];
        stat_GF1av [d][l] += stat_GF1 [d][l][count];
        stat_GF2av [d][l] += stat_GF2 [d][l][count];

      }

      stat_aGF1av [d][l] = stat_aGF1av [d][l]/(Nsamples[l]*BMsamples);
      stat_aGF2av [d][l] = stat_aGF2av [d][l]/(Nsamples[l]*BMsamples);
      stat_hybGFav [d][l] = stat_hybGFav [d][l]/(Nsamples[l]*BMsamples);
      stat_GF1av [d][l] = stat_GF1av [d][l]/(Nsamples[l]*BMsamples);
      stat_GF2av [d][l] = stat_GF2av [d][l]/(Nsamples[l]*BMsamples);
    }
    for ( int count=0; count<BMsamples; count++){
      stat_BMav [l] += stat_BM [l][count];
    }
    stat_BMav[l] = stat_BMav[l]/BMsamples;
  }

  double arrTav[6][nL];
  for (int l=0; l<nL; l++){
    for ( int d=0; d<7; d++){
      arrTav[d][l]=0;
    }
  }
  for ( int l=0; l<nL; l++ ){
    for ( int count=0; count<Nsamples[l]*BMsamples; count++){
      arrTav[0][l] += arrT [0][l][count];
      arrTav[1][l] += arrT [1][l][count];
      arrTav[2][l] += arrT [2][l][count];
      arrTav[3][l] += arrT [3][l][count];
      arrTav[4][l] += arrT [4][l][count];
    }
    arrTav [0][l] = arrTav [0][l]/(Nsamples[l]*BMsamples);
    arrTav [1][l] = arrTav [1][l]/(Nsamples[l]*BMsamples);
    arrTav [2][l] = arrTav [2][l]/(Nsamples[l]*BMsamples);
    arrTav [3][l] = arrTav [3][l]/(Nsamples[l]*BMsamples);
    arrTav [4][l] = arrTav [4][l]/(Nsamples[l]*BMsamples);

    for ( int count=0; count<BMsamples; count++){
      arrTav[5][l] += arrT [5][l][count];
    }
    arrTav [5][l] = arrTav [5][l]/BMsamples;
 
  }

  // standard deviations computing
  double stat_aGF1sd[3][nL], stat_aGF2sd[3][nL], stat_hybGFsd[3][nL], stat_GF1sd[3][nL],stat_GF2sd[3][nL],stat_BMsd[nL];
  for ( int l=0; l<nL; l++) {
    for ( int d=0; d<3; d++ ){
      stat_GF1sd[d][l] = 0;
      stat_GF2sd[d][l] = 0;
      stat_aGF1sd [d][l] = 0;
      stat_aGF2sd [d][l] = 0;
      stat_hybGFsd [d][l] = 0;
    }
    stat_BMsd [l] =0;
  }
  for ( int l=0; l<nL; l++) {
    for ( int d=0; d<3; d++){
      for ( int count=0; count<Nsamples[l]*BMsamples; count++){
        stat_aGF1sd [d][l] += pow (stat_aGF1av[d][l] - stat_aGF1 [d][l][count],2);
        stat_aGF2sd [d][l] += pow (stat_aGF2av[d][l] - stat_aGF2 [d][l][count],2);
        stat_hybGFsd [d][l] += pow (stat_hybGFav[d][l] - stat_hybGF [d][l][count],2);
        stat_GF1sd [d][l] += pow ( stat_GF1av[d][l] - stat_GF1 [d][l][count],2);
        stat_GF2sd [d][l] += pow ( stat_GF2av[d][l] - stat_GF2 [d][l][count],2);
      }
      stat_aGF1sd [d][l] = sqrt(stat_aGF1sd [d][l]/pow(Nsamples[l]*BMsamples,2));
      stat_aGF2sd [d][l] = sqrt(stat_aGF2sd [d][l]/pow(Nsamples[l]*BMsamples,2));
      stat_hybGFsd [d][l] = sqrt(stat_hybGFsd [d][l]/pow(Nsamples[l]*BMsamples,2));
      stat_GF1sd [d][l] = sqrt(stat_GF1sd [d][l]/pow(Nsamples[l]*BMsamples,2));
      stat_GF2sd [d][l] = sqrt(stat_GF2sd [d][l]/pow(Nsamples[l]*BMsamples,2));
    }
    for ( int count=0; count<BMsamples; count++){
      stat_BMsd [l] += pow( stat_BMav[l] - stat_BM [l][count],2);
    }
    stat_BMsd[l] = sqrt(stat_BMsd[l]/BMsamples/BMsamples);
  }

  double arrTsd[6][nL];
  for (int l=0; l<nL; l++){
    for ( int d=0; d<7; d++){
      arrTsd[d][l]=0;
    }
  }
  for ( int l=0; l<nL; l++ ){
    for ( int count=0; count<Nsamples[l]*BMsamples; count++){
      arrTsd[0][l] += pow (arrT [0][l][count] - arrTav[0][l],2);
      arrTsd[1][l] += pow (arrT [1][l][count] - arrTav[1][l],2);
      arrTsd[2][l] += pow (arrT [2][l][count] - arrTav[2][l],2);
      arrTsd[3][l] += pow (arrT [3][l][count] - arrTav[3][l],2);
      arrTsd[4][l] += pow (arrT [4][l][count] - arrTav[4][l],2);
    }      
    arrTsd [0][l] = sqrt(arrTsd [0][l]/pow(Nsamples[l]*BMsamples,2));
    arrTsd [1][l] = sqrt(arrTsd [1][l]/pow(Nsamples[l]*BMsamples,2));
    arrTsd [2][l] = sqrt(arrTsd [2][l]/pow(Nsamples[l]*BMsamples,2));
    arrTsd [3][l] = sqrt(arrTsd [3][l]/pow(Nsamples[l]*BMsamples,2));
    arrTsd [4][l] = sqrt(arrTsd [4][l]/pow(Nsamples[l]*BMsamples,2));

    for ( int count=0; count<BMsamples; count++){
      arrTsd[5][l] += pow (arrT [5][l][count] - arrTav[5][l],2);
    }
    arrTsd [5][l] = sqrt (arrTsd [5][l]/BMsamples/BMsamples);
 
  }


  for ( int l=0; l<nL; l++){

      printf( "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t",
              L[l], arrTav[0][l], arrTav[1][l], arrTav[2][l], arrTav[3][l], arrTav[4][l], arrTav[5][l],
              arrTsd [0][l], arrTsd[1][l], arrTsd[2][l], arrTsd[3][l], arrTsd[4][l], arrTsd[5][l] );

  }
    printf("\n\n");


  for ( int l=0; l<nL; l++){
      printf("%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t",
             L[l], stat_aGF1av[0][l],stat_aGF2av[0][l], stat_hybGFav[0][l], stat_GF1av[0][l], stat_GF2av[0][l],
             stat_aGF1sd[0][l], stat_aGF2sd[0][l],stat_hybGFsd[0][l],stat_GF1sd[0][l],stat_GF2sd[0][l]);
  }
    printf("\n\n");

    for ( int l=0; l<nL; l++){
        printf("%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t",
               L[l], stat_aGF1av[1][l],stat_aGF2av[1][l], stat_hybGFav[1][l], stat_GF1av[1][l], stat_GF2av[1][l],
               stat_aGF1sd[1][l], stat_aGF2sd[1][l],stat_hybGFsd[1][l],stat_GF1sd[1][l],stat_GF2sd[1][l]);
    }
    printf("\n\n");

    for ( int l=0; l<nL; l++){
        printf("%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t",
               L[l], stat_aGF1av[2][l],stat_aGF2av[2][l], stat_hybGFav[2][l], stat_GF1av[2][l], stat_GF2av[2][l],
               stat_aGF1sd[2][l], stat_aGF2sd[2][l],stat_hybGFsd[2][l],stat_GF1sd[2][l],stat_GF2sd[2][l]);
    }
    printf("\n\n");





}
