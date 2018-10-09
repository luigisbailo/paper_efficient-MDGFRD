//TO COMPILE: g++ -std=c++11 main.cpp -o main -lgsl -lgslcblas -lm

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_math.h>
#include <cmath>
#include <gsl/gsl_errno.h>
#include <iostream>
#include <fstream>
#include <algorithm> 
#include <vector>
#include <iomanip>
#include <chrono>
#include <cstring>
#include <sstream>

using namespace std::chrono;


#include "../parameters.hpp"
#include "../tools.hpp"
#include "../greensFunct.hpp"
#include "../draw.hpp"
#include "../init.hpp"
#include "../step.hpp"
#include "../shell.hpp"
#include "../print.hpp"
#include "../burst.hpp"
#include "../bruteForce.hpp"
#include "../checks.hpp"
#include "../run_GF.hpp"
#include "../run_aGF1.hpp"
#include "../run_aGF2.hpp"
#include "../run_hybGF.hpp"
#include "../run_BM.hpp"

int main ( int argc, char *argv[] ) {         

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
    Nsamples[0]=1000;
    Nsamples[0]=1000;
    Nsamples[0]=100;
    Nsamples[0]=100;
    Nsamples[0]=10;
    Nsamples[0]=1;

  double D_A=0.01;
  double D_B=0.01;
  double R_A=2.5;
  double R_B=2.5;
  double tau_bm=0.1;
  long int nINTsteps;
  int POWsteps=7;
  const int N = 10; 
  const int N_A = 5;
  const int N_B = 5;
  double alpha = 9.;
  int BMsamples = 2;



  high_resolution_clock::time_point t2,t1;
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
        for (int n=0; n<N; n++)
         diffStat[n];

        t1 = high_resolution_clock::now();

        run_aGF1 ( N_A, N_B, R_A, R_B, D_A, D_B, tau_bm, alpha, Tsim, L[l], stat, diffStat );

        t2 = high_resolution_clock::now();
        
        t12 = double(std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count())/1000000;

        arrT [0][l][count1*Nsamples[l]+count2] = t12;

        for (int d=0; d<3; d++)
          stat_aGF1 [d][l][count1*Nsamples[l]+count2] = stat[d]; 



        for (int d=0; d<3; d++)
         stat[d] = 0;
        for (int n=0; n<N; n++)
         diffStat[n];

        t1 = high_resolution_clock::now();

        run_aGF2 ( N_A, N_B, R_A, R_B, D_A, D_B, tau_bm, alpha, Tsim, L[l], stat, diffStat );

        t2 = high_resolution_clock::now();
        
        t12 = double(std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count())/1000000;

        arrT [1][l][count1*Nsamples[l]+count2] = t12;

        for (int d=0; d<3; d++)
          stat_aGF2 [d][l][count1*Nsamples[l]+count2] = stat[d]; 



        for (int d=0; d<3; d++)
         stat[d] = 0;
        for (int n=0; n<N; n++)
         diffStat[n];

        t1 = high_resolution_clock::now();

        run_hybGF ( N_A, N_B, R_A, R_B, D_A, D_B, tau_bm, alpha, Tsim, L[l], stat, diffStat );

        t2 = high_resolution_clock::now();
        
        t12 = double(std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count())/1000000;

        arrT [2][l][count1*Nsamples[l]+count2] = t12;

        for (int d=0; d<3; d++)
          stat_hybGF [d][l][count1*Nsamples[l]+count2] = stat[d]; 


        for (int d=0; d<3; d++)
          stat[d] = 0;
        for (int n=0; n<N; n++)
          diffStat[n];

        t1 = high_resolution_clock::now();

        run_GF ( N_A, N_B, R_A, R_B, D_A, D_B, 1.,1., tau_bm, alpha, Tsim, L[l], stat, diffStat );

        t2 = high_resolution_clock::now();
        
        t12 = double(std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count())/1000000;


        arrT [3][l][count1*Nsamples[l]+count2] = t12;

        for (int d=0; d<3; d++)
          stat_GF1 [d][l][count1*Nsamples[l]+count2] = stat[d]; 
      


        for (int d=0; d<3; d++)
          stat[d] = 0;
        for (int n=0; n<N; n++)
          diffStat[n];

        t1 = high_resolution_clock::now();

        run_GF ( N_A, N_B, R_A, R_B, D_A, D_B, 1.5, 2.5, tau_bm, alpha, Tsim, L[l], stat, diffStat );

        t2 = high_resolution_clock::now();
        
        t12 = double(std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count())/1000000;


        arrT [4][l][count1*Nsamples[l]+count2] = t12;

        for (int d=0; d<3; d++)
          stat_GF2 [d][l][count1*Nsamples[l]+count2] = stat[d]; 
      }



      for (int n=0; n<N; n++)
         diffStat[n];

      t1 = high_resolution_clock::now();

      run_BM ( N_A, N_B, R_A, R_B, D_A, D_B, tau_bm, Tsim, L[l], diffStat );

      t2 = high_resolution_clock::now();
      
 
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

        // cout << stat_aGF [d][l][count] << "--";
        // cout << stat_GF1 [d][l][count] << "\t"; 

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
      // cout << arrT[0][l][count] << endl;
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





  // // averages computing
  // double stat_aGFav[3][nL], stat_GFav[3][nL],stat_BMav[nL];
  // for ( int l=0; l<nL; l++) {
  //   for ( int d=0; d<3; d++ ){
  //     stat_GFav[d][l] = 0;
  //     stat_aGFav [d][l] = 0;
  //   }
  //   stat_BMav [l] =0;
  // }
  // for ( int l=0; l<nL; l++) {
  //   for ( int d=0; d<3; d++){

  //     for ( int count=0; count<Nsamples[l]*BMsamples; count++){
  //       stat_aGFav [d][l] += stat_aGF [d][l][count];
  //       stat_GFav [d][l] += stat_GF [d][l][count];

  //       // cout << stat_aGF [d][l][count] << "--";
  //       // cout << stat_GF [d][l][count] << "\t"; 

  //     }

  //     stat_aGFav [d][l] = stat_aGFav [d][l]/(Nsamples[l]*BMsamples);
  //     stat_GFav [d][l] = stat_GFav [d][l]/(Nsamples[l]*BMsamples);
  //   }
  //   for ( int count=0; count<BMsamples; count++){
  //     stat_BMav [l] += stat_BM [l][count];
  //   }
  //   stat_BMav[l] = stat_BMav[l]/BMsamples;
  // }

  // double arrTav[3][nL];
  // for (int l=0; l<nL; l++){
  //   for ( int d=0; d<4; d++){
  //     arrTav[d][l]=0;
  //   }
  // }
  // for ( int l=0; l<nL; l++ ){
  //   for ( int count=0; count<Nsamples[l]*BMsamples; count++){
  //     arrTav[0][l] += arrT [0][l][count];
  //     arrTav[1][l] += arrT [1][l][count];
  //   }      
  //   arrTav [0][l] = arrTav [0][l]/(Nsamples[l]*BMsamples);
  //   arrTav [1][l] = arrTav [1][l]/(Nsamples[l]*BMsamples);

  //   for ( int count=0; count<BMsamples; count++){
  //     arrTav[2][l] += arrT [2][l][count];
  //   }
  //   arrTav [2][l] = arrTav [2][l]/BMsamples;
 
  // }

  // // standard deviations computing
  // double stat_aGFsd[3][nL], stat_GFsd[3][nL],stat_BMsd[nL];
  // for ( int l=0; l<nL; l++) {
  //   for ( int d=0; d<3; d++ ){
  //     stat_GFsd[d][l] = 0;
  //     stat_aGFsd [d][l] = 0;
  //   }
  //   stat_BMsd [l] =0;
  // }
  // for ( int l=0; l<nL; l++) {
  //   for ( int d=0; d<3; d++){
  //     for ( int count=0; count<Nsamples[l]*BMsamples; count++){
  //       stat_aGFsd [d][l] += pow (stat_aGFav[d][l] - stat_aGF [d][l][count],2);
  //       stat_GFsd [d][l] += pow ( stat_GFav[d][l] - stat_GF [d][l][count],2);
  //     }
  //     stat_aGFsd [d][l] = sqrt(stat_aGFsd [d][l]/pow(Nsamples[l]*BMsamples,2));
  //     stat_GFsd [d][l] = sqrt(stat_GFsd [d][l]/pow(Nsamples[l]*BMsamples,2));
  //   }
  //   for ( int count=0; count<BMsamples; count++){
  //     stat_BMsd [l] += pow( stat_BMav[l] - stat_BM [l][count],2);
  //   }
  //   stat_BMsd[l] = sqrt(stat_BMsd[l]/BMsamples/BMsamples);
  // }

  // double arrTsd[3][nL];
  // for (int l=0; l<nL; l++){
  //   for ( int d=0; d<4; d++){
  //     arrTsd[d][l]=0;
  //   }
  // }
  // for ( int l=0; l<nL; l++ ){
  //   for ( int count=0; count<Nsamples[l]*BMsamples; count++){
  //     arrTsd[0][l] += pow (arrT [0][l][count] - arrTav[0][l],2);
  //     arrTsd[1][l] += pow (arrT [1][l][count] - arrTav[1][l],2);
  //   }      
  //   arrTsd [0][l] = sqrt(arrTsd [0][l]/pow(Nsamples[l]*BMsamples,2));
  //   arrTsd [1][l] = sqrt(arrTsd [1][l]/pow(Nsamples[l]*BMsamples,2));

  //   for ( int count=0; count<BMsamples; count++){
  //     arrTsd[2][l] += pow (arrT [2][l][count] - arrTav[2][l],2);
  //   }
  //   arrTsd [2][l] = sqrt (arrTsd [2][l]/BMsamples/BMsamples);
 
  // }



  // cout << setprecision(5)<< fixed; 

  // for ( int l=0; l<nL; l++){

  //   cout << int(L[l]) << "\t" << arrTav [0][l] << "\t" << arrTav[1][l] << "\t" << arrTav[2][l];    
  //   cout << "\t" << arrTsd [0][l] << "\t" << arrTsd[1][l] << "\t" << arrTsd[2][l] << endl;

  // }
  // cout << endl;

  // cout << setprecision(1); 

  // for ( int l=0; l<nL; l++){
  //   cout << int(L[l]) << "\t" << stat_aGFav[0][l] << "\t" << stat_GFav[0][l];
  //   cout << "\t" << stat_aGFsd[0][l] << "\t" << stat_GFsd[0][l] << endl;
  // }
  // cout << endl;

  // for ( int l=0; l<nL; l++){
  //   cout << int(L[l]) << "\t" << stat_aGFav[1][l] << "\t" << stat_GFav[1][l];
  //   cout << "\t" << stat_aGFsd [1][l] << "\t" << stat_GFsd[1][l] << endl;
  // }
  // cout << endl;


  // for ( int l=0; l<nL; l++){
  //   cout << int(L[l]) << "\t" << stat_aGFav[2][l] << "\t" << stat_GFav[2][l];
  //   cout << "\t" << stat_aGFsd [2][l] << "\t" << stat_GFsd[2][l] << endl;
  // }
  // cout << endl;

  // // for ( int l=0; l<nL; l++){
  // //   cout << int(L[l]) << "\t" << stat_aGFav[3][l] << "\t" << stat_GFav [3][l] << "\t" << stat_BMav[l] ;
  // //   cout << "\t" << stat_aGFsd[3][l] << "\t" << stat_GFsd [3][l] << "\t" << stat_BMsd[l]  << endl;
  // // }


 

}
