//TO COMPILE: g++ -std=c++11 main.cpp -o main -lgsl -lgslcblas -lm
using namespace std;
#define print(x) cout << x << endl;

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

using namespace std;
using namespace std::chrono;


#include "../../2draft/parameters.hpp"
#include "../tools.hpp"
#include "../../2draft/greensFunct.hpp"
#include "../../2draft/draw.hpp"
#include "../init.hpp"
#include "../step.hpp"
#include "../shell.hpp"
#include "../print.hpp"
#include "../burst.hpp"
#include "../bruteForce.hpp"
#include "../checks.hpp"
#include "../run_GF.hpp"
#include "../run_aGF.hpp"
#include "../run_BM.hpp"

int main ( int argc, char *argv[] ) {         

  // int nL = 4;
  // int Nsamples[nL]={100,10,10,1};
  // int L[nL] = {5,50,500,5000}; 


  int nN = 2;
  int Nsamples[nN]={5,1};
  // int Nsamples[nN]={1,1,1,1};
  int N[nN] = {1000,10000}; 

  double D_A=0.01;
  double D_B=0.01;
  double R_A=2.5;
  double R_B=2.5;
  double tau_bm=0.1;
  double alpha = 9.;
  int BMsamples = 1;

  long int nINTsteps;
  double L;
  int POWsteps;
  double maxSh_aGF[nN];
  double maxSh_GF[nN];
  double BMgrid;

  stringstream convert_L (argv[1]); 
  stringstream convert_POWSTEPS (argv[2]); 


  if (!(convert_L >> L ))
    exit (EXIT_FAILURE);  
  if (!(convert_POWSTEPS >> POWsteps ))
    exit (EXIT_FAILURE);  

//  if (L==200){
    maxSh_GF[0]=5; 
    maxSh_GF[1]=5; 
    maxSh_GF[2]=5; 
    maxSh_GF[3]=5; 
    maxSh_aGF[0]=2.5; 
    maxSh_aGF[1]=2.5; 
    maxSh_aGF[2]=2.5; 
    maxSh_aGF[3]=2.5;
    BMgrid = 5; 
//  }
//  else if (L=2100){
//    maxSh_GF[0]=10; 
//    maxSh_GF[1]=10; 
//    maxSh_GF[2]=10; 
//    maxSh_GF[3]=10; 
//    maxSh_aGF[0]=10; 
//    maxSh_aGF[1]=10; 
//    maxSh_aGF[2]=10; 
//    maxSh_aGF[3]=10;
//    BMgrid = 25; 
//  }
//  else {
//    cout << "ERROR: box size L";
//    exit(EXIT_FAILURE);
//  }

  high_resolution_clock::time_point t2,t1;
  double t12, avT;

  nINTsteps = pow (10,POWsteps);
  const double Tsim = nINTsteps*tau_bm;

  int stat [3];
  // double diffStat[N];

  // stat[0] number burstings
  // stat[1] number domains
  // stat[2] BD steps
  // stat[3] squared displacement
  double stat_aGF[3][nN][Nsamples[0]*BMsamples];
  double stat_GF1[3][nN][Nsamples[0]*BMsamples];
  double stat_GF2[3][nN][Nsamples[0]*BMsamples];
  double stat_BM [nN][BMsamples];
  double arrT [4][nN][Nsamples[0]*BMsamples];

  cout << L << endl;
  cout << D_A << endl;
  cout << D_B << endl;
  cout << tau_bm << endl;
  cout << nINTsteps << endl;
  cout << alpha << endl;
  cout << "Minimal aGF shell: " << alpha*sqrt(D_A*tau_bm) << "\t" << alpha*sqrt(D_B*tau_bm) << endl;
  cout << "Minimal GF1 shell:  " << 1*R_A << "\t" << 1*R_B << endl;
  cout << "Minimal GF2 shell:  " << 1.5*R_A << "\t" << 2.5*R_B << endl;
  cout << endl;


  for ( int count1=0; count1 < BMsamples; count1 ++){

    int n=0;

    for ( int l=0; l<nN; l++) {

      for ( int count2=0; count2<Nsamples[l]; count2++ ) {

        for (int d=0; d<3; d++)
         stat[d] = 0;
        // for (int n=0; n<N; n++)
        //  diffStat[n];

        t1 = high_resolution_clock::now();

        run_aGF ( N[l], 0, R_A, R_B, D_A, D_B, tau_bm, alpha, Tsim, L, maxSh_aGF[l], stat );

        t2 = high_resolution_clock::now();
        
        t12 = double(std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count())/1000000;

        arrT [0][l][count1*Nsamples[l]+count2] = t12;

        for (int d=0; d<3; d++)
          stat_aGF [d][l][count1*Nsamples[l]+count2] = stat[d]; 


        for (int d=0; d<3; d++)
          stat[d] = 0;
        // for (int n=0; n<N; n++)
        //   diffStat[n];

        t1 = high_resolution_clock::now();

        // if ( L>300 | N[l]<5000){
          // run_GF ( N[l], 0, R_A, R_B, D_A, D_B, 1.,1., tau_bm, alpha, Tsim, L, maxSh_GF[l], stat);
        // }
        
        t2 = high_resolution_clock::now();
        
        t12 = double(std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count())/1000000;


        arrT [1][l][count1*Nsamples[l]+count2] = t12;

        for (int d=0; d<3; d++)
          stat_GF1 [d][l][count1*Nsamples[l]+count2] = stat[d]; 
      

      }



      // for (int n=0; n<N; n++)
      //    diffStat[n];


      t1 = high_resolution_clock::now();

      run_BM ( N[l], 0, R_A, R_B, D_A, D_B, tau_bm, Tsim, L, BMgrid);

      t2 = high_resolution_clock::now();
      
      t12 = double(std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count())/1000000;

      arrT [3][l][count1] = t12;

      stat_BM [l][count1] = stat[3]; 

      }


  }



  // averages computing
  double stat_aGFav[3][nN], stat_GF1av[3][nN],stat_GF2av[3][nN],stat_BMav[nN];
  for ( int l=0; l<nN; l++) {
    for ( int d=0; d<3; d++ ){
      stat_GF1av[d][l] = 0;
      stat_GF2av[d][l] = 0;
      stat_aGFav [d][l] = 0;
    }
    stat_BMav [l] =0;
  }
  for ( int l=0; l<nN; l++) {
    for ( int d=0; d<3; d++){

      for ( int count=0; count<Nsamples[l]*BMsamples; count++){
        stat_aGFav [d][l] += stat_aGF [d][l][count];
        stat_GF1av [d][l] += stat_GF1 [d][l][count];
        // stat_GF2av [d][l] += stat_GF2 [d][l][count];

        // cout << stat_aGF [d][l][count] << "--";
        // cout << stat_GF [d][l][count] << "\t"; 

      }

      stat_aGFav [d][l] = stat_aGFav [d][l]/(Nsamples[l]*BMsamples);
      stat_GF1av [d][l] = stat_GF1av [d][l]/(Nsamples[l]*BMsamples);
      // stat_GF2av [d][l] = stat_GF2av [d][l]/(Nsamples[l]*BMsamples);
    }
    for ( int count=0; count<BMsamples; count++){
      stat_BMav [l] += stat_BM [l][count];
    }
    stat_BMav[l] = stat_BMav[l]/BMsamples;
  }

  double arrTav[3][nN];
  for (int l=0; l<nN; l++){
    for ( int d=0; d<4; d++){
      arrTav[d][l]=0;
    }
  }
  for ( int l=0; l<nN; l++ ){
    for ( int count=0; count<Nsamples[l]*BMsamples; count++){
      arrTav[0][l] += arrT [0][l][count];
      arrTav[1][l] += arrT [1][l][count];
      // arrTav[2][l] += arrT [2][l][count];
    }      
    arrTav [0][l] = arrTav [0][l]/(Nsamples[l]*BMsamples);
    arrTav [1][l] = arrTav [1][l]/(Nsamples[l]*BMsamples);
    // arrTav [2][l] = arrTav [2][l]/(Nsamples[l]*BMsamples);

    for ( int count=0; count<BMsamples; count++){
      arrTav[3][l] += arrT [3][l][count];
    }
    arrTav [3][l] = arrTav [3][l]/BMsamples;
 
  }

  // standard deviations computing
  double stat_aGFsd[3][nN], stat_GF1sd[3][nN],stat_GF2sd[3][nN],stat_BMsd[nN];
  for ( int l=0; l<nN; l++) {
    for ( int d=0; d<3; d++ ){
      stat_GF1sd[d][l] = 0;
      // stat_GF2sd[d][l] = 0;
      stat_aGFsd [d][l] = 0;
    }
    stat_BMsd [l] =0;
  }
  for ( int l=0; l<nN; l++) {
    for ( int d=0; d<3; d++){
      for ( int count=0; count<Nsamples[l]*BMsamples; count++){
        stat_aGFsd [d][l] += pow (stat_aGFav[d][l] - stat_aGF [d][l][count],2);
        stat_GF1sd [d][l] += pow ( stat_GF1av[d][l] - stat_GF1 [d][l][count],2);
        // stat_GF2sd [d][l] += pow ( stat_GF2av[d][l] - stat_GF2 [d][l][count],2);
      }
      stat_aGFsd [d][l] = sqrt(stat_aGFsd [d][l]/pow(Nsamples[l]*BMsamples,2));
      stat_GF1sd [d][l] = sqrt(stat_GF1sd [d][l]/pow(Nsamples[l]*BMsamples,2));
      // stat_GF2sd [d][l] = sqrt(stat_GF2sd [d][l]/pow(Nsamples[l]*BMsamples,2));
    }
    for ( int count=0; count<BMsamples; count++){
      stat_BMsd [l] += pow( stat_BMav[l] - stat_BM [l][count],2);
    }
    stat_BMsd[l] = sqrt(stat_BMsd[l]/BMsamples/BMsamples);
  }

  double arrTsd[3][nN];
  for (int l=0; l<nN; l++){
    for ( int d=0; d<4; d++){
      arrTsd[d][l]=0;
    }
  }
  for ( int l=0; l<nN; l++ ){
    for ( int count=0; count<Nsamples[l]*BMsamples; count++){
      arrTsd[0][l] += pow (arrT [0][l][count] - arrTav[0][l],2);
      arrTsd[1][l] += pow (arrT [1][l][count] - arrTav[1][l],2);
      // arrTsd[2][l] += pow (arrT [2][l][count] - arrTav[2][l],2);
    }      
    arrTsd [0][l] = sqrt(arrTsd [0][l]/pow(Nsamples[l]*BMsamples,2));
    arrTsd [1][l] = sqrt(arrTsd [1][l]/pow(Nsamples[l]*BMsamples,2));
    // arrTsd [2][l] = sqrt(arrTsd [2][l]/pow(Nsamples[l]*BMsamples,2));

    for ( int count=0; count<BMsamples; count++){
      arrTsd[3][l] += pow (arrT [3][l][count] - arrTav[3][l],2);
    }
    arrTsd [3][l] = sqrt (arrTsd [3][l]/BMsamples/BMsamples);
 
  }



  cout << setprecision(5)<< fixed; 

  for ( int l=0; l<nN; l++){

    cout << int(N[l]) << "\t" << arrTav [0][l] << "\t" << arrTav[1][l] << "\t" << arrTav[3][l];    
    cout << "\t" << arrTsd [0][l] << "\t" << arrTsd[1][l]  << "\t" << arrTsd[3][l] << endl;

  }
  cout << endl;

  cout << setprecision(1); 

  for ( int l=0; l<nN; l++){
    cout << int(N[l]) << "\t" << stat_aGFav[0][l] << "\t" << stat_GF1av[0][l] ;
    cout << "\t" << stat_aGFsd[0][l] << "\t" << stat_GF1sd[0][l]  << endl;
  }
  cout << endl;

  for ( int l=0; l<nN; l++){
    cout << int(N[l]) << "\t" << stat_aGFav[1][l] << "\t" << stat_GF1av[1][l];
    cout << "\t" << stat_aGFsd [1][l] << "\t" << stat_GF1sd[1][l]  << endl;
  }
  cout << endl;


  for ( int l=0; l<nN; l++){
    cout << int(N[l]) << "\t" << stat_aGFav[2][l] << "\t" << stat_GF1av[2][l] ;
    cout << "\t" << stat_aGFsd [2][l] << "\t" << stat_GF1sd[2][l] << endl;
  }
  cout << endl;





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
