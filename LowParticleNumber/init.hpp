// author luigisbailo


void initPos_aGF1 ( particle *particles, gsl_rng *r, int N_A, int N_B, double R_A, double R_B, double D_A, double D_B, double tau_bm, double alpha, double L) {

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


    if (count < N_A){        

      particles[count].type = 0;
      particles[count].Diff = D_A;    
      particles[count].sqrtDiff = sqrt(D_A);       
      particles[count].radius = R_A;
      particles[count].R_bd = alpha*sqrt(D_A*tau_bm); 
      particles[count].shellRed = 0;
      particles[count].burstR = particles[count].R_bd; 

    }
    else{

      particles[count].type = 1;
      particles[count].Diff = D_B;    
      particles[count].sqrtDiff = sqrt(D_B); 
      particles[count].radius = R_B;
      particles[count].R_bd = alpha*sqrt(D_B*tau_bm); 
      particles[count].shellRed = 0;
      particles[count].burstR = particles[count].R_bd; 

     }
  }
}


void initPos_aGF2 ( particle *particles, gsl_rng *r, int N_A, int N_B, double R_A, double R_B, double D_A, double D_B, double tau_bm, double alpha, double L) {

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


    if (count < N_A){        

      particles[count].type = 0;
      particles[count].Diff = D_A;    
      particles[count].sqrtDiff = sqrt(D_A);       
      particles[count].radius = R_A;
      particles[count].R_bd = alpha*sqrt(D_A*tau_bm); 
      particles[count].shellRed = 5*sqrt(2*D_A*tau_bm);
      particles[count].burstR = particles[count].R_bd; 

    }
    else{

      particles[count].type = 1;
      particles[count].Diff = D_B;    
      particles[count].sqrtDiff = sqrt(D_B); 
      particles[count].radius = R_B;
      particles[count].R_bd = alpha*sqrt(D_B*tau_bm); 
      particles[count].shellRed = 5*sqrt(2*D_B*tau_bm);
      particles[count].burstR = particles[count].R_bd; 

     }
  }
}


void initPos_hybGF ( particle *particles, gsl_rng *r, int N_A, int N_B, double R_A, double R_B, double D_A, double D_B, double tau_bm, double alpha, double L) {

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


    if (count < N_A){        

      particles[count].type = 0;
      particles[count].Diff = D_A;    
      particles[count].sqrtDiff = sqrt(D_A); 
      particles[count].radius = R_A;
      particles[count].R_gfrd = alpha*sqrt(D_A*tau_bm);
      particles[count].R_bd = alpha*sqrt(D_A*tau_bm);
      particles[count].burstR = alpha*sqrt(std::min(D_A,D_B)*tau_bm);
      
    }
    else{

      particles[count].type = 1;
      particles[count].Diff = D_B;    
      particles[count].sqrtDiff = sqrt(D_B); 
      particles[count].radius = R_B;
      particles[count].R_gfrd = alpha*sqrt(D_B*tau_bm); 
      particles[count].R_bd = alpha*sqrt(D_A*tau_bm); 
      particles[count].burstR = alpha*sqrt(std::min(D_A,D_B)*tau_bm);

    }
  }
}


void initPos_GF ( particle *particles, gsl_rng *r, int N_A, int N_B, double R_A, double R_B, double D_A, double D_B, double MU_BM, double MU_GF, double tau_bm, double alpha, double L) {

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


    if (count < N_A){        

      particles[count].type = 0;
      particles[count].Diff = D_A;    
      particles[count].sqrtDiff = sqrt(D_A); 
      particles[count].radius = R_A;
      particles[count].R_gfrd = MU_GF*R_A;
      particles[count].R_bd = MU_BM*R_A;
      particles[count].burstR = MU_BM*std::min(R_A,R_B);
      
    }
    else{

      particles[count].type = 1;
      particles[count].Diff = D_B;    
      particles[count].sqrtDiff = sqrt(D_B); 
      particles[count].radius = R_B;
      particles[count].R_gfrd = MU_GF*R_B; 
      particles[count].R_bd = MU_BM*R_B; 
      particles[count].burstR = MU_BM*std::min(R_A,R_B);

    }
  }
}


void initPos_BM ( particle *particles, gsl_rng *r, int N_A, int N_B, double R_A, double R_B, double D_A, double D_B, double tau_bm, double L) {

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



void initShell_aGF ( particle *particles, gsl_rng *r, int N, double tau_bm, double sqrt2TAU_BM, double L, int *stat ) {

  double initShells [N];

  for (int i=0; i<N; i++){

    for (int j=0; j<N; j++) {

      if (i!=j){
        initShells[j] = (sqrt (dist2_per (&particles[i],&particles[j],L) ) - particles[i].radius - particles[j].radius ) /
              ( 1 +  particles[j].sqrtDiff / particles[i].sqrtDiff  );    
      }
      else {
        initShells[j] = L;
      }

    }
    
  double R = *std::min_element( initShells, initShells+N);


  if ( R > particles[i].R_bd ){

    (*stat) ++;
    double deltaPos [3];

    polarTransf ( deltaPos, R, gsl_rng_uniform (r), gsl_rng_uniform (r));

    particles[i].tau_exit = drawTimeNewt ( R, particles[i].Diff, gsl_rng_uniform(r) );
    particles[i].shell = R;
    particles[i].gf = true;
    particles[i].pos_exit[0]= particles[i].pos[0] + deltaPos[0];
    particles[i].pos_exit[1]= particles[i].pos[1] + deltaPos[1];
    particles[i].pos_exit[2]= particles[i].pos[2] + deltaPos[2];
    checkBound ( particles[i].pos_exit, particles[i].pos_period, L );

  }      
  else {
    
    particles[i].tau_exit += tau_bm;
    particles[i].pos_exit[0] += gsl_ran_gaussian (r,1) * particles[i].sqrtDiff * sqrt2TAU_BM;
    particles[i].pos_exit[1] += gsl_ran_gaussian (r,1) * particles[i].sqrtDiff * sqrt2TAU_BM;
    particles[i].pos_exit[2] += gsl_ran_gaussian (r,1) * particles[i].sqrtDiff * sqrt2TAU_BM;
    checkBound ( particles[i].pos_exit, particles[i].pos_period, L );

  }

 }

}


void initShell_GF ( particle *particles, gsl_rng *r, int N, double tau_bm, double sqrt2TAU_BM, double L, int *stat ) {

  double initShells [N];

  for (int i=0; i<N; i++){

    for (int j=0; j<N; j++) {

      if (i!=j){
        initShells[j] = (sqrt (dist2_per (&particles[i],&particles[j],L) ) - particles[i].radius - particles[j].radius ) / 2;
      }
      else {
        initShells[j] = L;
      }

    }
    
  double R = *std::min_element( initShells, initShells+N);


  if ( R > particles[i].R_bd ){

    (*stat) ++;
 
   particles[i].tau_exit = drawTimeNewt ( R, particles[i].Diff, gsl_rng_uniform(r) ); 
    particles[i].shell = R;
    particles[i].gf = true;
    particles[i].pos_exit[0]= -1;
    particles[i].pos_exit[1]= -1;
    particles[i].pos_exit[2]= -1;

  }      
  else {

    particles[i].tau_exit += tau_bm;
    particles[i].pos_exit[0] += gsl_ran_gaussian (r,1)*particles[i].sqrtDiff * sqrt2TAU_BM;
    particles[i].pos_exit[1] += gsl_ran_gaussian (r,1)*particles[i].sqrtDiff * sqrt2TAU_BM;
    particles[i].pos_exit[2] += gsl_ran_gaussian (r,1)*particles[i].sqrtDiff * sqrt2TAU_BM;
    checkBound (particles[i].pos_exit,particles[i].pos_period, L);

  }

 }

}

