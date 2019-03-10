// author luigisbailo


void fig5_8 ( double D_A, double D_B, double R_A, double R_B ) {

    int nL = 7;
    double L[nL];
    L[0]=2500;
    L[1]=1150;
    L[2]=545;
    L[3]=255;
    L[4]=118;
    L[5]=55;
    L[6]=25;
    int Nsamples[nL];
    Nsamples[0]=1000;
    Nsamples[1]=1000;
    Nsamples[2]=1000;
    Nsamples[3]=100;
    Nsamples[4]=100;
    Nsamples[5]=10;
    Nsamples[6]=1;

  double tau_bm=0.1;
  long int nINTsteps;
  int POWsteps=5;
  const int N = 10; 
  const int N_A = 5;
  const int N_B = 5;
  double alpha = 9.;
  int BMsamples = 2;



  std::chrono::high_resolution_clock::time_point t2,t1;
  double t12, avT;

  nINTsteps = int(pow (10,POWsteps));
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

  std::cout << D_A << std::endl;
    std::cout << D_B << std::endl;
    std::cout << tau_bm << std::endl;
    std::cout << nINTsteps << std::endl;
    std::cout << alpha << std::endl;
    std::cout << "Expected squared distance per particle: " << 6*Tsim*(N_A*D_A+N_B*D_B)/N << std::endl;
    std::cout << "Minimal aGF shell: " << alpha*sqrt(D_A*tau_bm) << "\t" << alpha*sqrt(D_B*tau_bm) << std::endl;
    std::cout << "Minimal GF1 shell:  " << 1*R_A << "\t" << 1*R_B << std::endl;
    std::cout << "Minimal GF2 shell:  " << 1.5*R_A << "\t" << 2.5*R_B << std::endl;
    std::cout << std::endl;


  for ( int count1=0; count1 < BMsamples; count1 ++){

    int n=0;

    for ( int l=0; l<nL; l++) {

      for ( int count2=0; count2<Nsamples[l]; count2++ ) {


        for (int d=0; d<3; d++)
         stat[d] = 0;

        t1 = std::chrono::high_resolution_clock::now();

        run_aGF1 ( N_A, N_B, R_A, R_B, D_A, D_B, tau_bm, alpha, Tsim, L[l], stat, diffStat );

        t2 = std::chrono::high_resolution_clock::now();
        
        t12 = double(std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count())/1000000;

        arrT [0][l][count1*Nsamples[l]+count2] = t12;

        for (int d=0; d<3; d++)
          stat_aGF1 [d][l][count1*Nsamples[l]+count2] = stat[d]; 



        for (int d=0; d<3; d++)
         stat[d] = 0;
        for (int n=0; n<N; n++)
         diffStat[n]= 0;

        t1 = std::chrono::high_resolution_clock::now();

        run_aGF2 ( N_A, N_B, R_A, R_B, D_A, D_B, tau_bm, alpha, Tsim, L[l], stat, diffStat );

        t2 = std::chrono::high_resolution_clock::now();
        
        t12 = double(std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count())/1000000;

        arrT [1][l][count1*Nsamples[l]+count2] = t12;

        for (int d=0; d<3; d++)
          stat_aGF2 [d][l][count1*Nsamples[l]+count2] = stat[d]; 



        for (int d=0; d<3; d++)
         stat[d] = 0;
        for (int n=0; n<N; n++)
         diffStat[n]=0;

        t1 = std::chrono::high_resolution_clock::now();

        run_hybGF ( N_A, N_B, R_A, R_B, D_A, D_B, tau_bm, alpha, Tsim, L[l], stat, diffStat );

        t2 = std::chrono::high_resolution_clock::now();
        
        t12 = double(std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count())/1000000;

        arrT [2][l][count1*Nsamples[l]+count2] = t12;

        for (int d=0; d<3; d++)
          stat_hybGF [d][l][count1*Nsamples[l]+count2] = stat[d]; 


        for (int d=0; d<3; d++)
          stat[d] = 0;
        for (int n=0; n<N; n++)
          diffStat[n] = 0;

        t1 = std::chrono::high_resolution_clock::now();

        run_GF ( N_A, N_B, R_A, R_B, D_A, D_B, 1.,1., tau_bm, alpha, Tsim, L[l], stat, diffStat );

        t2 = std::chrono::high_resolution_clock::now();
        
        t12 = double(std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count())/1000000;


        arrT [3][l][count1*Nsamples[l]+count2] = t12;

        for (int d=0; d<3; d++)
          stat_GF1 [d][l][count1*Nsamples[l]+count2] = stat[d]; 
      


        for (int d=0; d<3; d++)
          stat[d] = 0;
        for (int n=0; n<N; n++)
          diffStat[n] = 0;

        t1 = std::chrono::high_resolution_clock::now();

        run_GF ( N_A, N_B, R_A, R_B, D_A, D_B, 1.5, 2.5, tau_bm, alpha, Tsim, L[l], stat, diffStat );

        t2 = std::chrono::high_resolution_clock::now();
        
        t12 = double(std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count())/1000000;


        arrT [4][l][count1*Nsamples[l]+count2] = t12;

        for (int d=0; d<3; d++)
          stat_GF2 [d][l][count1*Nsamples[l]+count2] = stat[d]; 
      }



      for (int n=0; n<N; n++)
         diffStat[n] = 0;

      t1 = std::chrono::high_resolution_clock::now();

      run_BM ( N_A, N_B, R_A, R_B, D_A, D_B, tau_bm, Tsim, L[l], diffStat );

      t2 = std::chrono::high_resolution_clock::now();
      
 
      arrT [5][l][count1] = t12;

      stat_BM [l][count1] = stat[3]; 

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



    std::cout << std::setprecision(5)<< std::fixed;

  for ( int l=0; l<nL; l++){

      std::cout << int(L[l]) << "\t" << arrTav [0][l] << "\t" << arrTav[1][l] << "\t" << arrTav[2][l] << "\t" << arrTav[3][l] << "\t" << arrTav[4][l] << "\t" << arrTav[5][l];
      std::cout << "\t" << arrTsd [0][l] << "\t" << arrTsd[1][l] << "\t" << arrTsd[2][l] << "\t" << arrTsd[3][l] << "\t" << arrTsd[4][l] << "\t" << arrTsd[5][l] << std::endl;

  }
    std::cout <<std:: endl;

    std::cout << std::setprecision(1);

  for ( int l=0; l<nL; l++){
      std::cout << int(L[l]) << "\t" << stat_aGF1av[0][l] << "\t" << stat_aGF2av[0][l] << "\t" << stat_hybGFav[0][l] << "\t" << stat_GF1av[0][l] << "\t" << stat_GF2av[0][l];
      std::cout << "\t" << stat_aGF1sd[0][l] << "\t" << stat_aGF2sd[0][l] << "\t" <<stat_hybGFsd[0][l] << "\t" << stat_GF1sd[0][l] << "\t" << stat_GF2sd[0][l] << std::endl;
  }
    std::cout << std::endl;

  for ( int l=0; l<nL; l++){
      std::cout << int(L[l]) << "\t" << stat_aGF1av[1][l]<< "\t" << stat_aGF2av[1][l] << "\t" << stat_hybGFav[1][l] << "\t" << stat_GF1av[1][l]<< "\t" << stat_GF2av[1][l];
      std::cout << "\t" << stat_aGF1sd [1][l] << "\t" << stat_aGF2sd [1][l] << "\t" << stat_hybGFsd [1][l] << "\t" << stat_GF1sd[1][l] << "\t" << stat_GF2sd[1][l] << std::endl;
  }
    std::cout << std::endl;


  for ( int l=0; l<nL; l++){
      std::cout << int(L[l]) << "\t" << stat_aGF1av[2][l] << "\t" << stat_aGF2av[2][l] << "\t" << stat_hybGFav[2][l] << "\t" << stat_GF1av[2][l] << "\t" << stat_GF2av[2][l];
      std::cout << "\t" << stat_aGF1sd [2][l] << "\t" << stat_aGF2sd [2][l] << "\t"  << stat_hybGFsd [2][l] << "\t" << stat_GF1sd[2][l] << "\t" << stat_GF2sd[2][l] << std::endl;
  }
    std::cout << std::endl;


}
