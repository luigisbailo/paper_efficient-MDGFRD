// author luigisbailo


void initPos_aGF_nl ( struct particle_nl *particles, gsl_rng *r, int N_A, int N_B, double R_A, double R_B, double D_A, double D_B,
                      double tau_bm, double alpha, double L, double boxSizeGF ) {

  double x,y,z;

  for (int count=0; count < N_A+N_B; count++){

    particles[count].label = count;
    particles[count].time = 0;
    particles[count].tau_exit = 0; 
    particles[count].shell = 0; 
    particles[count].burst = false; 
    particles[count].gf = false;

    x = gsl_rng_uniform (r)*L;
    y = gsl_rng_uniform (r)*L;  
    z = gsl_rng_uniform (r)*L;

    particles[count].pos[0] = x;
    particles[count].pos[1] = y;
    particles[count].pos[2] = z;   
    particles[count].pos_exit[0] = x;
    particles[count].pos_exit[1] = y;
    particles[count].pos_exit[2] = z; 
    particles[count].pos_init[0] = x;
    particles[count].pos_init[1] = y;
    particles[count].pos_init[2] = z; 
    particles[count].pos_period[0] = 0;
    particles[count].pos_period[1] = 0;
    particles[count].pos_period[2] = 0;
    particles[count].idxBox [0] = (int) trunc(x/boxSizeGF);
    particles[count].idxBox [1] = (int) trunc(y/boxSizeGF);
    particles[count].idxBox [2] = (int) trunc(z/boxSizeGF);

    if (count < N_A){        

      particles[count].type = 0;
      particles[count].Diff = D_A;    
      particles[count].sqrtDiff = sqrt(D_A);       
      particles[count].radius = R_A;
      particles[count].R_bd = alpha*sqrt(D_A*tau_bm); 
      particles[count].shellRed = 4*sqrt(2*D_A*tau_bm);
      particles[count].burstR = particles[count].R_bd; 

    }
    else{

      particles[count].type = 1;
      particles[count].Diff = D_B;    
      particles[count].sqrtDiff = sqrt(D_B); 
      particles[count].radius = R_B;
      particles[count].R_bd = alpha*sqrt(D_B*tau_bm); 
      particles[count].shellRed = 4*sqrt(2*D_A*tau_bm);
      particles[count].burstR = particles[count].R_bd; 

     }
  }
}


void initPos_GF_nl ( struct particle_nl *particles, gsl_rng *r, int N_A, int N_B, double R_A, double R_B, double D_A, double D_B,
                     double MU_BM, double MU_GF, double tau_bm, double alpha, double L, double boxSizeGF ) {

  double x,y,z;

  for (int count=0; count < N_A+N_B; count++){

    particles[count].label = count;
    particles[count].time = 0;
    particles[count].tau_exit = 0; 
    particles[count].shell = 0; 
    particles[count].burst = false; 
    particles[count].gf = false;

    x = gsl_rng_uniform (r)*L;
    y = gsl_rng_uniform (r)*L;  
    z = gsl_rng_uniform (r)*L;

    particles[count].pos[0] = x;
    particles[count].pos[1] = y;
    particles[count].pos[2] = z;   
    particles[count].pos_exit[0] = x;
    particles[count].pos_exit[1] = y;
    particles[count].pos_exit[2] = z;   
    particles[count].pos_init[0] = x;
    particles[count].pos_init[1] = y;
    particles[count].pos_init[2] = z; 
    particles[count].pos_period[0] = 0;
    particles[count].pos_period[1] = 0;
    particles[count].pos_period[2] = 0;
    particles[count].idxBox [0] = trunc(x/boxSizeGF);
    particles[count].idxBox [1] = trunc(y/boxSizeGF);
    particles[count].idxBox [2] = trunc(z/boxSizeGF);

    if (count < N_A){        

      particles[count].type = 0;
      particles[count].Diff = D_A;    
      particles[count].sqrtDiff = sqrt(D_A); 
      particles[count].radius = R_A;
      // particles[count].R_gfrd = alpha*sqrt(D_A*tau_bm); 
      // particles[count].R_bd = alpha*sqrt(D_A*tau_bm); 
      particles[count].R_gfrd = MU_GF*R_A;
      particles[count].R_bd = MU_BM*R_A;
      // particles[count].burstR = particles[count].R_bd; 
      particles[count].burstR = MU_BM*fmin(R_A,R_B);
      
    }
    else{

      particles[count].type = 1;
      particles[count].Diff = D_B;    
      particles[count].sqrtDiff = sqrt(D_B); 
      particles[count].radius = R_B;
      // particles[count].R_gfrd = alpha*sqrt(D_B*tau_bm); 
      // particles[count].R_bd = alpha*sqrt(D_B*tau_bm); 
      particles[count].R_gfrd = MU_GF*R_B; 
      particles[count].R_bd = MU_BM*R_B; 
      // particles[count].burstR = particles[count].R_bd;    
      particles[count].burstR = MU_BM*fmin(R_A,R_B);

    }
  }
}


void initPos_BM_nl ( struct particle_nl *particles, gsl_rng *r, int N_A, int N_B, double R_A, double R_B, double D_A, double D_B,
                     double tau_bm, double L, double boxSize) {

  double x,y,z;

  for (int count=0; count < N_A+N_B; count++){

    particles[count].label = count;
    particles[count].time = 0;
    particles[count].tau_exit = 0; 
    particles[count].shell = 0; 
    particles[count].burst = false; 
    particles[count].gf = false;

    x = gsl_rng_uniform (r)*L;
    y = gsl_rng_uniform (r)*L;  
    z = gsl_rng_uniform (r)*L;

    particles[count].pos[0] = x;
    particles[count].pos[1] = y;
    particles[count].pos[2] = z;   
    particles[count].pos_exit[0] = x;
    particles[count].pos_exit[1] = y;
    particles[count].pos_exit[2] = z;   
    particles[count].pos_init[0] = x;
    particles[count].pos_init[1] = y;
    particles[count].pos_init[2] = z; 
    particles[count].pos_period[0] = 0;
    particles[count].pos_period[1] = 0;
    particles[count].pos_period[2] = 0;
    particles[count].idxBox [0] = trunc(x/boxSize);
    particles[count].idxBox [1] = trunc(y/boxSize);
    particles[count].idxBox [2] = trunc(z/boxSize);


    if (count < N_A){        
      particles[count].Diff = D_A;    
      particles[count].sqrtDiff = sqrt(D_A); 
      particles[count].radius = R_A;
      
    }
    else{

      particles[count].Diff = D_B;    
      particles[count].sqrtDiff = sqrt(D_B); 
      particles[count].radius = R_B;


    }
  }
}



void initShell_aGF_nl ( struct particle_nl *particles, struct boxcell *grid, double *distRow, int *distLabel,
                        gsl_rng *r, int N, double tau_bm, double sqrt2TAU_BM, double maxShell, double L, int nBoxes,
                        int *stat ) {


  double Shell,R,Rmin,Dij[2];
  int j;

  for (int iPart=0; iPart<N; iPart++){

    int ilabel = iPart;
    getDist_nl (particles, grid, distRow, distLabel, ilabel, N, L, nBoxes );

    Rmin = particles[iPart].R_bd;
    Shell = maxShell;
    Dij[0] = particles[iPart].sqrtDiff;
    j=0;

    // it loops over all particles to check whether they are within the bursting radius
    while (!(distLabel[j]<0)){

      int jPart = distLabel[j];
      Dij[1] = particles[jPart].sqrtDiff;

      R = distRow[j] / ( 1 +  Dij[1] / Dij[0]  );

      if (R<Shell)
        Shell = R ;
      if (Shell<Rmin){
        Shell = 0;
        break;
      }
      j++;

    }
    if ( Shell > particles[iPart].R_bd ){

      (*stat) ++;
      double deltaPos [3];

      polarTransf_nl ( deltaPos, Shell, gsl_rng_uniform (r), gsl_rng_uniform (r));

      particles[iPart].tau_exit = drawTimeNewt ( Shell, particles[iPart].Diff, gsl_rng_uniform(r) ); 
      particles[iPart].shell = Shell;
      particles[iPart].gf = true;
      particles[iPart].pos_exit[0]= particles[iPart].pos[0] + deltaPos[0];
      particles[iPart].pos_exit[1]= particles[iPart].pos[1] + deltaPos[1];
      particles[iPart].pos_exit[2]= particles[iPart].pos[2] + deltaPos[2];
      checkBound_nl ( particles[iPart].pos_exit, particles[iPart].pos_period, L );

    }
    else {
      
      particles[iPart].tau_exit += tau_bm;
      particles[iPart].pos_exit[0] += gsl_ran_gaussian (r,1) * particles[iPart].sqrtDiff * sqrt2TAU_BM;
      particles[iPart].pos_exit[1] += gsl_ran_gaussian (r,1) * particles[iPart].sqrtDiff * sqrt2TAU_BM;
      particles[iPart].pos_exit[2] += gsl_ran_gaussian (r,1) * particles[iPart].sqrtDiff * sqrt2TAU_BM;
      checkBound_nl ( particles[iPart].pos_exit, particles[iPart].pos_period, L );

    }

  }

}


void initShell_GF_nl ( struct particle_nl *particles, struct boxcell *grid, double *distRow, int *distLabel,
                       gsl_rng *r, int N, double tau_bm, double sqrt2TAU_BM, double maxShell, double L, int nBoxes,
                       int *stat ) {


  double Shell,R,Rmin,Dij[2];
  int j;

  for (int iPart=0; iPart<N; iPart++){

    int ilabel = iPart;
    getDist_nl (particles, grid, distRow, distLabel, ilabel, N, L, nBoxes );

    Rmin = particles[iPart].R_bd;
    Shell = maxShell;
    Dij[0] = particles[iPart].sqrtDiff;
    j=0;

    // it cycles over all particles to check whether they are within the bursting radius
    while (!(distLabel[j]<0)){

      int jPart = distLabel[j];
      Dij[1] = particles[jPart].sqrtDiff;

      R = distRow[j] / 2;

      if (R<Shell)
        Shell = R ;
      if (Shell<Rmin){
        Shell = 0;
        break;
      }
      j++;

    }
    

    if ( Shell > particles[iPart].R_bd ){

      (*stat) ++;
      double deltaPos [3];

      polarTransf_nl ( deltaPos, Shell, gsl_rng_uniform (r), gsl_rng_uniform (r));

      particles[iPart].tau_exit = drawTimeNewt ( Shell, particles[iPart].Diff, gsl_rng_uniform(r) ); 
      particles[iPart].shell = Shell;
      particles[iPart].gf = true;
      particles[iPart].pos_exit[0]= particles[iPart].pos[0] + deltaPos[0];
      particles[iPart].pos_exit[1]= particles[iPart].pos[1] + deltaPos[1];
      particles[iPart].pos_exit[2]= particles[iPart].pos[2] + deltaPos[2];
      checkBound_nl ( particles[iPart].pos_exit, particles[iPart].pos_period, L );

    }      
    else {
      
      particles[iPart].tau_exit += tau_bm;
      particles[iPart].pos_exit[0] += gsl_ran_gaussian (r,1) * particles[iPart].sqrtDiff * sqrt2TAU_BM;
      particles[iPart].pos_exit[1] += gsl_ran_gaussian (r,1) * particles[iPart].sqrtDiff * sqrt2TAU_BM;
      particles[iPart].pos_exit[2] += gsl_ran_gaussian (r,1) * particles[iPart].sqrtDiff * sqrt2TAU_BM;
      checkBound_nl ( particles[iPart].pos_exit, particles[iPart].pos_period, L );

    }

  }

}



void initGridGF_nl (struct particle_nl *particles, struct boxcell *grid,
                    int nBoxesGF, double boxSizeGF, int N ) {

  for ( int i=0; i<nBoxesGF; i++){
    for ( int j=0; j<nBoxesGF; j++){
      for ( int k=0; k<nBoxesGF; k++){

        grid [i*nBoxesGF*nBoxesGF + j*nBoxesGF +k].element[0]=0;
        for ( int n=1; n<NLIMBOX; n++){

          grid [i*nBoxesGF*nBoxesGF + j*nBoxesGF +k].element[n]=-1;

        }          
      }
    }
  }


  int i,j,k;
  for ( int n=0; n<N; n++){

    i=(int)trunc(particles[n].pos[0]/boxSizeGF);
    j=(int)trunc(particles[n].pos[1]/boxSizeGF);
    k=(int)trunc(particles[n].pos[2]/boxSizeGF);
    grid[i*nBoxesGF*nBoxesGF + j*nBoxesGF +k].element[0]++;
    grid[i*nBoxesGF*nBoxesGF + j*nBoxesGF +k].element[grid[i*nBoxesGF*nBoxesGF + j*nBoxesGF +k].element[0]]=n;

  }

}


void initGridBM_nl ( struct particle_nl *particles, struct boxcell *grid,
                     int nBoxes, double boxSize, int N){


  for ( int i=0; i<nBoxes; i++){
    for ( int j=0; j<nBoxes; j++){
      for ( int k=0; k<nBoxes; k++){
        grid [i*nBoxes*nBoxes + j*nBoxes +k].element[0]=0;
        for ( int n=1; n<NLIMBOX; n++){
          grid [i*nBoxes*nBoxes + j*nBoxes +k].element[n]=-1;
        }          

      }
    }
  }

  int i,j,k;
  for ( int n=0; n<N; n++){

    i=trunc(particles[n].pos[0]/boxSize);
    j=trunc(particles[n].pos[1]/boxSize);
    k=trunc(particles[n].pos[2]/boxSize);

    grid[i*nBoxes*nBoxes + j*nBoxes +k].element[0]++;
    grid[i*nBoxes*nBoxes + j*nBoxes +k].element[grid[i*nBoxes*nBoxes + j*nBoxes +k].element[0]]=n;

  }

}
