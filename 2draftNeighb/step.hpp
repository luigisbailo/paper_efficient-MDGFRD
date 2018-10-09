void GFstep_aGF ( particle *myPart, gsl_rng *r, double R, double L ){
//The shell radius "R" has been already determined with getR () and it has been already checked that the domain can be constructed

  double deltaPos [3];

  //extraction of the exit time and exit position in polar coordinates
  polarTransf ( deltaPos, R, gsl_rng_uniform (r), gsl_rng_uniform (r) );

  //deltaPos now contains the displacements in cartesian coordinates
  myPart->pos_exit[0] += deltaPos[0];  
  myPart->pos_exit[1] += deltaPos[1];
  myPart->pos_exit[2] += deltaPos[2];      
  checkBound (myPart->pos_exit,myPart->pos_period, L );
  myPart->tau_exit += drawTimeNewt ( R, myPart->Diff, gsl_rng_uniform(r) );
  myPart->shell = R;

    
}


void GFstep_GF ( particle *myPart, gsl_rng *r, double R ){
//The shell radius "R" has been already determined with getR () and it has been already checked that the domain can be constructed

  //extraction of the exit time and exit position in polar coordinates
  myPart->tau_exit += drawTimeNewt ( R, myPart->Diff, gsl_rng_uniform(r) );
  myPart->pos_exit[0] = -1;
  myPart->pos_exit[1] = -1;
  myPart->pos_exit[2] = -1;
  myPart->shell = R;

    
}



void BMstep ( particle *particles, int *partList, double *distRow, int *distLabel, gsl_rng *r, double tau_bm, double sqrt2TAU_BM, int N, double L ) {
//sqrtTAU_BM is sqrt(2*TAU_BM)

  double dist,deltaPos[3], varPos[3];
  int j = 0;

  deltaPos[0] = 0;
  deltaPos[1] = 0;
  deltaPos[2] = 0;

  while ( !(distLabel[j]<0) ){


    int jPart = distLabel[j];


    dist = distRow[j];

    if ( dist < 0 ){
      // cout << "INT with part: " << distLabel[j] << endl;
      // if the distance with another particle is lower than R_INTER we take into account their interaction
      // varPos is the cartesian projection of the particles distancce 
      // the origin is centered in the count position

      varPos[0] = particles[partList[0]].pos[0] - particles[jPart].pos[0];
        if (varPos[0]>L/2) varPos[0] -= L;
        else if (varPos[0]<-L/2) varPos[0] += L;
      varPos[1] = particles[partList[0]].pos[1] - particles[jPart].pos[1];
        if (varPos[1]>L/2) varPos[1] -= L;
        else if (varPos[1]<-L/2) varPos[1] += L;
      varPos[2] = particles[partList[0]].pos[2] - particles[jPart].pos[2];
        if (varPos[2]>L/2) varPos[2] -= L;
        else if (varPos[2]<-L/2) varPos[2] += L;



      // The particles interact via a repulsive harmonic interaction
      deltaPos[0] += K*particles[partList[0]].Diff * (varPos[0]/(dist+particles[partList[0]].radius+particles[jPart].radius)) * (-dist) * tau_bm;
      deltaPos[1] += K*particles[partList[0]].Diff * (varPos[1]/(dist+particles[partList[0]].radius+particles[jPart].radius)) * (-dist) * tau_bm;
      deltaPos[2] += K*particles[partList[0]].Diff * (varPos[2]/(dist+particles[partList[0]].radius+particles[jPart].radius)) * (-dist) * tau_bm;

    }

    j++;  

  }

  particles[partList[0]].pos_exit[0] = deltaPos[0] + particles[partList[0]].pos[0] + gsl_ran_gaussian (r,1)* particles[partList[0]].sqrtDiff * sqrt2TAU_BM;
  particles[partList[0]].pos_exit[1] = deltaPos[1] + particles[partList[0]].pos[1] + gsl_ran_gaussian (r,1)* particles[partList[0]].sqrtDiff * sqrt2TAU_BM;    
  particles[partList[0]].pos_exit[2] = deltaPos[2] + particles[partList[0]].pos[2] + gsl_ran_gaussian (r,1)* particles[partList[0]].sqrtDiff * sqrt2TAU_BM;
  checkBound (particles[partList[0]].pos_exit, particles[partList[0]].pos_period, L );
  particles[partList[0]].tau_exit += tau_bm;
     
}



void synchPart_P_aGF ( particle *particles, int *partList, gsl_rng *r, int N, double Tsynch, double L ) {

  double synchPos [3];

  for ( int n=0; n<N; n++) {

    if (particles[n].shell>0) {

      
      if ( Tsynch < particles[n].time ){

         cout << "ERROR: synch" << "\n";
         // exit (EXIT_FAILURE);
         cout << setprecision(5);
         printPos_per ( particles, partList, N );

      }


        if (Tsynch-particles[n].time> (particles[n].shell*particles[n].shell)/particles[n].Diff/100){


          double Rsynch =  drawPosNewt ( Tsynch -particles[n].time, particles[n].shell, particles[n].Diff, gsl_rng_uniform(r) );
          polarTransf ( synchPos, Rsynch, gsl_rng_uniform (r), gsl_rng_uniform (r) );

          fixBound_aGF ( particles[n].pos, particles[n].pos_exit, particles[n].pos_period, L);        
          particles[n].pos[0] += synchPos[0];
          particles[n].pos[1] += synchPos[1];
          particles[n].pos[2] += synchPos[2];
          checkBound (particles[n].pos,particles[n].pos_period, L );
          particles[n].shell = Rsynch;
          particles[n].time = Tsynch;
      
        }

        else {

          double sqrt2dt = sqrt (2*(Tsynch-particles[n].time));

          fixBound_aGF ( particles[n].pos, particles[n].pos_exit, particles[n].pos_period, L);        
          particles[n].pos[0] += gsl_ran_gaussian (r,1)*particles[n].sqrtDiff * sqrt2dt;
          particles[n].pos[1] += gsl_ran_gaussian (r,1)*particles[n].sqrtDiff * sqrt2dt;
          particles[n].pos[2] += gsl_ran_gaussian (r,1)*particles[n].sqrtDiff * sqrt2dt;
          checkBound (particles[n].pos,particles[n].pos_period, L );
          particles[n].shell = sqrt2dt;
          particles[n].time = Tsynch;
 

        }

    }
  
  }
    
}




void synchPart_P_GF ( particle *particles, int *partList, gsl_rng *r, int N, double Tsynch, double L ) {

  double synchPos [3];

  for ( int n=0; n<N; n++) {

    if (particles[n].shell>0) {

      
      if ( Tsynch < particles[n].time ){

         cout << "ERROR: synch" << "\n";
         cout << setprecision(5);
    printPos_per ( particles, partList, N );
         // exit (EXIT_FAILURE);

      }


        if (Tsynch-particles[n].time> (particles[n].shell*particles[n].shell)/particles[n].Diff/100){


          double Rsynch =  drawPosNewt ( Tsynch -particles[n].time, particles[n].shell, particles[n].Diff, gsl_rng_uniform(r) );
          polarTransf ( synchPos, Rsynch, gsl_rng_uniform (r), gsl_rng_uniform (r));

          particles[n].pos[0] += synchPos[0];
          particles[n].pos[1] += synchPos[1];
          particles[n].pos[2] += synchPos[2];
          checkBound (particles[n].pos,particles[n].pos_period, L );
          particles[n].shell = Rsynch;
          particles[n].time = Tsynch;
      
        }

        else {

          double sqrt2dt = sqrt (2*(Tsynch-particles[n].time));
          particles[n].pos[0] += gsl_ran_gaussian (r,1)*particles[n].sqrtDiff * sqrt2dt;
          particles[n].pos[1] += gsl_ran_gaussian (r,1)*particles[n].sqrtDiff * sqrt2dt;
          particles[n].pos[2] += gsl_ran_gaussian (r,1)*particles[n].sqrtDiff * sqrt2dt;
          checkBound (particles[n].pos,particles[n].pos_period, L );
          particles[n].shell = sqrt2dt;
          particles[n].time = Tsynch;
 

        }

    }
  
  }
    
}
