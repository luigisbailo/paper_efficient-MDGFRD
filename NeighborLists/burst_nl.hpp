// author luigisbailo


void burst_P_aGF_nl ( struct particle_nl *particles, int *partList, double *distRow, int *distLabel,
                      gsl_rng *r,  int N, double L ) {

  double deltaPos [3]; 
  int iPart = partList[0];
  int j=0;

  // it cycles over all particles to check whether they are within the bursting radius
  while (!(distLabel[j]<0)){


    int jPart = distLabel[j];

          
    if (  particles[jPart].gf == true  &&  distRow[j] - particles[jPart].shell < particles[iPart].burstR ) {

      particles[jPart].burst = true;  
      particles[jPart].gf = false;


      //The P function is not sampled at very small times, when the survival function S can be approximated to 1      
      if (particles[iPart].time-particles[jPart].time >
              (particles[jPart].shell*particles[jPart].shell)/particles[jPart].Diff/1000){

        polarTransf_nl ( deltaPos, drawPosNewt ( particles[iPart].time-particles[jPart].time,
                                                 particles[jPart].shell, particles[jPart].Diff, gsl_rng_uniform(r) ),
                     gsl_rng_uniform (r), gsl_rng_uniform (r) );
        //deltaPos now contains the displacements in cartesian coordinates

        fixBound_aGF_nl ( particles[jPart].pos, particles[jPart].pos_exit, particles[jPart].pos_period, L);
        particles[jPart].pos[0] += deltaPos[0];  
        particles[jPart].pos[1] += deltaPos[1];
        particles[jPart].pos[2] += deltaPos[2];      
        checkBound_nl ( particles[jPart].pos, particles[jPart].pos_period, L );
        particles[jPart].pos_exit[0]=particles[jPart].pos[0];
        particles[jPart].pos_exit[1]=particles[jPart].pos[1];
        particles[jPart].pos_exit[2]=particles[jPart].pos[2];
        particles[jPart].shell = 0;
        particles[jPart].time = particles[iPart].time;
        particles[jPart].tau_exit = particles[iPart].time;

      }
      else if (particles[iPart].time>particles[jPart].time) {

        //At very small times, the bursting procedure consists simply in a brownian motion integration step 
        double sqrt2dt = sqrt (2*(particles[iPart].time-particles[jPart].time));
        fixBound_aGF_nl ( particles[jPart].pos, particles[jPart].pos_exit, particles[jPart].pos_period, L);
        particles[jPart].pos[0] += gsl_ran_gaussian (r,1)*particles[jPart].sqrtDiff * sqrt2dt;
        particles[jPart].pos[1] += gsl_ran_gaussian (r,1)*particles[jPart].sqrtDiff * sqrt2dt;
        particles[jPart].pos[2] += gsl_ran_gaussian (r,1)*particles[jPart].sqrtDiff * sqrt2dt;
        checkBound_nl ( particles[jPart].pos, particles[jPart].pos_period, L );
        particles[jPart].pos_exit[0]=particles[jPart].pos[0];
        particles[jPart].pos_exit[1]=particles[jPart].pos[1];
        particles[jPart].pos_exit[2]=particles[jPart].pos[2];
        particles[jPart].shell = 0;
        particles[jPart].time = particles[iPart].time;
        particles[jPart].tau_exit = particles[iPart].time;


      }
      else{
        
        //In case the domain is burst at the same time of the construction
        fixBound_aGF_nl ( particles[jPart].pos, particles[jPart].pos_exit, particles[jPart].pos_period, L);
        particles[jPart].pos_exit[0]=particles[jPart].pos[0];
        particles[jPart].pos_exit[1]=particles[jPart].pos[1];
        particles[jPart].pos_exit[2]=particles[jPart].pos[2];
        particles[jPart].shell = 0;
        particles[jPart].tau_exit = particles[iPart].time;


      }

      // "distRow[]" is updated with the new distances,
      // and whether there is a new closest distance to insert in distRow[0] is checked
      distRow [j] = sqrt(dist2_per_nl ( &particles[iPart], &particles[jPart], L )) -
              particles[iPart].radius - particles[jPart].radius;
      sortBurst_nl (jPart,partList,N);

    }

    j ++;
  }  

}


void burst_P_GF_nl ( struct particle_nl *particles, int *partList, double *distRow, int* distLabel,
                     gsl_rng *r,  int N, double L ) {


  double deltaPos [3]; 
  int iPart = partList[0];
  int j=0;

  // it cycles over all particles to check whether they are within the bursting radius
  while (!(distLabel[j]<0)){


    int jPart = distLabel[j];

          
    if ( particles[jPart].gf == true && distRow[j] - particles[jPart].shell < particles[iPart].burstR ){

      particles[jPart].burst = true;  
      particles[jPart].gf = false;

      //The P function is not sampled at very small times, when the survival function S can be approximated to 1      
      if (particles[iPart].time-particles[jPart].time >
              (particles[jPart].shell*particles[jPart].shell)/particles[jPart].Diff/1000){
      
        polarTransf_nl ( deltaPos, drawPosNewt ( particles[iPart].time-particles[jPart].time,  particles[jPart].shell,
                                                 particles[jPart].Diff, gsl_rng_uniform(r) ),
                      gsl_rng_uniform (r), gsl_rng_uniform (r) );
        //deltaPos now contains the displacements in cartesian coordinates
        
        particles[jPart].pos[0] += deltaPos[0];  
        particles[jPart].pos[1] += deltaPos[1];
        particles[jPart].pos[2] += deltaPos[2];      
        checkBound_nl ( particles[jPart].pos, particles[jPart].pos_period, L );
        particles[jPart].pos_exit[0] = -1;
        particles[jPart].pos_exit[1] = -1;
        particles[jPart].pos_exit[2] = -1;
        particles[jPart].shell = 0;
        particles[jPart].time = particles[iPart].time;
        particles[jPart].tau_exit = particles[iPart].time;

      }
      else if (particles[iPart].time>particles[jPart].time){

        //At very small times, the bursting procedure consists simply in a brownian motion integration step 
        double sqrt2dt = sqrt (2*(particles[iPart].time-particles[jPart].time));
        particles[jPart].pos[0] += gsl_ran_gaussian (r,1)*particles[jPart].sqrtDiff * sqrt2dt;
        particles[jPart].pos[1] += gsl_ran_gaussian (r,1)*particles[jPart].sqrtDiff * sqrt2dt;
        particles[jPart].pos[2] += gsl_ran_gaussian (r,1)*particles[jPart].sqrtDiff * sqrt2dt;
        checkBound_nl ( particles[jPart].pos, particles[jPart].pos_period, L );
        particles[jPart].shell = 0;
        particles[jPart].time = particles[iPart].time;
        particles[jPart].tau_exit = particles[iPart].time;

      }
      else{
        
        //In case the domain is burst at the same time of the construction
        particles[jPart].pos_exit[0]=-1;
        particles[jPart].pos_exit[1]=-1;
        particles[jPart].pos_exit[2]=-1;
        particles[jPart].shell = 0;
        particles[jPart].tau_exit = particles[iPart].time;

      }

      // "distRow[]" is updated with the new distances, and weather there is a new closest distance to insert in distRow[0] is checked    
      distRow [j] = sqrt(dist2_per_nl ( &particles[iPart], &particles[jPart], L ))
                    - particles[iPart].radius - particles[jPart].radius;
      sortBurst_nl (jPart,partList,N);
    
    }

    j++;
  }  

}

