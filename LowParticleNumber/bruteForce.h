// author luigisbailo


void getBFdistances ( struct particle *particles, struct BFdistances *d, int N, double L){

    for ( int i=0; i<N; i++ ) {

        for ( int j=i+1; j<N; j++ ){

            d[i].dd[j] = sqrt(dist2_per(&particles[i],&particles[j],L));

        }
    }
}


void BFstep ( struct particle *particles, struct BFdistances *d, gsl_rng *r,
              double tau_bm, int N, double sqrt2TAU_BM, double L ) {
//dist,XYZ,deltaPos,varPos are just pointers to external free memory 
//sqrtTAU_BM is sqrt(2*TAU_BM)
//deltaPos is an array of the increments in the position
//particles are in a box modeled with a soft core repulsion on the boundaries

  double dist,deltaPos [3], varPos[3];

  for (int i=0; i<N; i++) {

    deltaPos[0] = 0;
    deltaPos[1] = 0;
    deltaPos[2] = 0;


    for( int j=0; j<N; j++){

      if (i==j) continue;
      // if the distance with another particle is lower than R_INTER we take into account their interaction
      if (i<j){
           dist = d[i].dd[j];
      }
      else{
           dist = d[j].dd[i];
      } 

      if ( dist <particles[i].radius+particles[j].radius){

        // varPos is the cartesian projection of the particles distancce 
        // the origin is centered in the count position
        varPos[0] = particles[i].pos[0]-particles[j].pos[0];
            if (varPos[0]>L/2) varPos[0] -= L;
            else if (varPos[0]<-L/2) varPos[0] += L;
        varPos[1] = particles[i].pos[1]-particles[j].pos[1];
            if (varPos[1]>L/2) varPos[1] -= L;
            else if (varPos[1]<-L/2) varPos[1] += L;
        varPos[2] = particles[i].pos[2]-particles[j].pos[2];
            if (varPos[2]>L/2) varPos[2] -= L;
            else if (varPos[2]<-L/2) varPos[2] += L;

        // The particles interact via a repulsive harmonic interaction
        deltaPos[0] += K*particles[i].Diff * (varPos[0]/ dist ) *
                (particles[i].radius+particles[j].radius -dist) * tau_bm;
        deltaPos[1] += K*particles[i].Diff * (varPos[1]/ dist ) *
                (particles[i].radius+particles[j].radius -dist) * tau_bm;
        deltaPos[2] += K*particles[i].Diff * (varPos[2]/ dist ) *
                (particles[i].radius+particles[j].radius -dist) * tau_bm;

        if ( fabs(sqrt(varPos[0]*varPos[0]+varPos[1]*varPos[1]+varPos[2]*varPos[2])- dist) > 0.001 ){
            printf ("ERROR: distances in forces computing in BM");
        }

      }

    }

    //brownian displacement
    deltaPos[0] += gsl_ran_gaussian (r,1)*particles[i].sqrtDiff * sqrt2TAU_BM;
    deltaPos[1] += gsl_ran_gaussian (r,1)*particles[i].sqrtDiff * sqrt2TAU_BM;
    deltaPos[2] += gsl_ran_gaussian (r,1)*particles[i].sqrtDiff * sqrt2TAU_BM;

    particles[i].pos_exit[0] += deltaPos[0];
    particles[i].pos_exit[1] += deltaPos[1];    
    particles[i].pos_exit[2] += deltaPos[2];
    particles[i].tau_exit += tau_bm;
    checkBound (particles[i].pos_exit, particles[i].pos_period, L ); 

  }     

}


void BFupdate ( struct particle *particles, int N) {

  for ( int n=0; n<N; n++){

    particles[n].pos[0] = particles[n].pos_exit[0];
    particles[n].pos[1] = particles[n].pos_exit[1];
    particles[n].pos[2] = particles[n].pos_exit[2];
    particles[n].time = particles[n].tau_exit;

  }

}