// author luigisbailo


double getR_aGF ( struct particle *particles, int *particleList, double *shells, double *distRow, int N, double L ){
//shell radius computed considering the expectd exit position of the other particles
//since it is considering only the exit positions, and not the actual shells it might overlap other shells
//this function is called after checking that the particle has space to construct a domain

  int iPart,jPart;
  double Dij[2];
  double Rshell[3];

  iPart = particleList[0];

  //Dij is an array of two elements, the first-[0] contains the diff coeff of our particle,
  // the second-[1] contains the diff coeff of neighbour particles
  Dij[0] = particles[iPart].Diff;

  //it cycles over all particles determining the shell radius predicted according to the aGF algorithm
  //the lowest radius is finally selected
  for ( int j=1; j<N; j++) {


    jPart = particleList [j];

    Dij[1] = particles[jPart].Diff;


    if ( particles[jPart].gf == false ) {


      if ( particles[jPart].type == particles[iPart].type )
        shells [j] = distRow [j] / 2;
      else
        shells [j] = distRow [j] / (1 + particles[jPart].sqrtDiff/particles[iPart].sqrtDiff );

      shells[j] = shells[j] - particles[jPart].shellRed;

    }
    else{

      Rshell[0] = distRow[j] - particles[jPart].shell;

      Rshell[2] = sqrt(dist2next_per(&particles[iPart],&particles[jPart],L))
                  - particles[iPart].radius - particles[jPart].radius;

      if ( particles[jPart].tau_exit - particles[iPart].time > Rshell[2]*Rshell[2]/6/Dij[0] ){  
      //if the first exit time of the j-part is higher than the average first exit time of the i-part from the largest domains,
      //the largest domain to i is proposed

        Rshell[1] = Rshell[2];

      }
      else if ( Dij[0] != Dij[1] ) { 

        double deltaRoot = (Rshell[2]*Rshell[2])/(36*Dij[1]*Dij[1]) - 
                          ((1/(6*Dij[1])-1/(6*Dij[0])) * ((Rshell[2]*Rshell[2]/(6*Dij[1]))
                                                        + particles[jPart].tau_exit - particles[iPart].time));
        Rshell[1] = ( Rshell[2]/(6*Dij[1]) - sqrt(deltaRoot) ) / (1/(6*Dij[1])-1/(6*Dij[0]));

      }      
      else {

        Rshell[1] = Rshell[2]/2 + (3 * Dij[0] * (particles[jPart].tau_exit - particles[iPart].time) ) / Rshell[2];       

      } 

      shells[j] = fmin (Rshell[0],Rshell[1]);

      if (particles[jPart].tau_exit-particles[iPart].time<shells[j]*shells[j]/6/Dij[0] ){
	    shells[j] = shells[j] - particles[jPart].shellRed;
	  }

    }

  }

  return min_element ( shells, 1, N );

}



double getR_GF ( struct particle *particles, int *particleList, double *shells, double *distRow, int N, double L ){

  int iPart,jPart;
  double Rshell;

  iPart = particleList[0];

  //shells[0] contains the minimial shell
  if ( particles[iPart].gf == true )  
    shells[0] = particles[iPart].R_bd; 
  else 
    shells[0] = particles[iPart].R_gfrd;

  for ( int j=1; j<N; j++) {

    jPart = particleList [j];

    // Rshell is an array of 3 elements,
    // the 1st-[0] contains the distance with the j-particle, the 2nd the distance with the expected exit position
    if ( particles[jPart].gf == false ){
      Rshell =  distRow [j] / 2;
    }
    else {
      Rshell = distRow [j] - particles[jPart].shell;
    }

    if( Rshell < shells[0] )
      return 0;
    else 
      shells[j] = Rshell;

  }

  return min_element( shells, 1, N);

}


