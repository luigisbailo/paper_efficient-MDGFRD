void BFstep ( particle *particles, vector<vector<int>> *grid, gsl_rng *r, int nBoxes, double boxSize, double tau_bm, int N, double sqrt2TAU_BM, double L ) {
//dist,XYZ,deltaPos,varPos are just pointers to external free memory 
//sqrtTAU_BM is sqrt(2*TAU_BM)
//deltaPos is an array of the increments in the position
//particles are in a box modeled with a soft core repulsion on the boundaries

  double dist,deltaPos [3], varPos[3];
  int iBoxi, iBoxj, iBoxk; 
  int iPart,jPart;
  double R;
  int i,j,k;

  for ( iPart=0; iPart<N; iPart++) {

    deltaPos[0] = 0;
    deltaPos[1] = 0;
    deltaPos[2] = 0;

    iBoxi = trunc (particles[iPart].pos[0]/boxSize);
    iBoxj = trunc (particles[iPart].pos[1]/boxSize);
    iBoxk = trunc (particles[iPart].pos[2]/boxSize);

    for (int ii=iBoxi-1; ii<iBoxi+2; ii++){
      i = ii;
      if (ii<0) i=nBoxes-1; 
      if (ii==nBoxes) i=0; 

      for (int jj=iBoxj-1; jj<iBoxj+2; jj++){
        j =jj;
        if (jj<0) j=nBoxes-1; 
        if (jj==nBoxes) j=0; 

        for (int kk=iBoxk-1; kk<iBoxk+2; kk++){
          k =kk;
          if (kk<0) k=nBoxes-1; 
          if (kk==nBoxes) k=0; 


          int nElem = (*grid)[i*nBoxes*nBoxes + j*nBoxes + k][0]; 
          int jElem=0;

          while (jElem<nElem){

            jElem++;
            jPart= (*grid)[i*nBoxes*nBoxes + j*nBoxes + k][jElem];
            if (iPart==jPart)
              continue;

            R = sqrt(dist2_per ( &particles[iPart], &particles[jPart], L ));

            if ( R <particles[iPart].radius+particles[jPart].radius){
              // varPos is the cartesian projection of the particles distancce c
              // the origin is centered in the count position
              varPos[0] = particles[iPart].pos[0]-particles[jPart].pos[0];
                  if (varPos[0]>L/2) varPos[0] -= L;
                  else if (varPos[0]<-L/2) varPos[0] += L;
              varPos[1] = particles[iPart].pos[1]-particles[jPart].pos[1];
                  if (varPos[1]>L/2) varPos[1] -= L;
                  else if (varPos[1]<-L/2) varPos[1] += L;
              varPos[2] = particles[iPart].pos[2]-particles[jPart].pos[2];
                  if (varPos[2]>L/2) varPos[2] -= L;
                  else if (varPos[2]<-L/2) varPos[2] += L;

              // The particles interact via a repulsive harmonic interaction
              deltaPos[0] += K*particles[iPart].Diff * (varPos[0]/ R ) * (particles[iPart].radius+particles[jPart].radius - R) * tau_bm;
              deltaPos[1] += K*particles[iPart].Diff * (varPos[1]/ R ) * (particles[iPart].radius+particles[jPart].radius - R) * tau_bm;
              deltaPos[2] += K*particles[iPart].Diff * (varPos[2]/ R ) * (particles[iPart].radius+particles[jPart].radius - R) * tau_bm;

            }
          }
        }
      }

    }

  
    //brownian displacement
    deltaPos[0] += gsl_ran_gaussian (r,1)*particles[iPart].sqrtDiff * sqrt2TAU_BM;
    deltaPos[1] += gsl_ran_gaussian (r,1)*particles[iPart].sqrtDiff * sqrt2TAU_BM;
    deltaPos[2] += gsl_ran_gaussian (r,1)*particles[iPart].sqrtDiff * sqrt2TAU_BM;


    particles[iPart].pos_exit[0] += deltaPos[0];
    particles[iPart].pos_exit[1] += deltaPos[1];    
    particles[iPart].pos_exit[2] += deltaPos[2];
    particles[iPart].tau_exit += tau_bm;
    checkBound (particles[iPart].pos_exit, particles[iPart].pos_period, L ); 

  }

}



void BFupdate ( particle *particles, int N) {


  for ( int n=0; n<N; n++){

    particles[n].pos[0] = particles[n].pos_exit[0];
    particles[n].pos[1] = particles[n].pos_exit[1];
    particles[n].pos[2] = particles[n].pos_exit[2];
    particles[n].time = particles[n].tau_exit;

  }

}


void updateGridBM ( particle *particles, vector <vector<int>> *grid, int nBoxes, double boxSize, int N){


  for ( int n=0; n<N; n++){

    int i,j,k;
    i = trunc (particles[n].pos[0]/boxSize);
    j = trunc (particles[n].pos[1]/boxSize);
    k = trunc (particles[n].pos[2]/boxSize);

// cout << i << "-" << j << "-" << k << "\t" << particles[n].pos[0] <<  endl;

    int nElem = (*grid) [i*nBoxes*nBoxes + j*nBoxes +k][0];
    int count = 0;
    while (count<nElem){
      count ++;
      (*grid) [i*nBoxes*nBoxes + j*nBoxes +k][count]=-1;
    }
    (*grid) [i*nBoxes*nBoxes + j*nBoxes +k][0]=0;

  }


  for ( int n=0; n<N; n++){

    int i,j,k;
    i=trunc(particles[n].pos[0]/boxSize);
    j=trunc(particles[n].pos[1]/boxSize);
    k=trunc(particles[n].pos[2]/boxSize);

    (*grid)[i*nBoxes*nBoxes + j*nBoxes +k][0]++;
    (*grid)[i*nBoxes*nBoxes + j*nBoxes +k][(*grid)[i*nBoxes*nBoxes + j*nBoxes +k][0]]=n;

  }

}
