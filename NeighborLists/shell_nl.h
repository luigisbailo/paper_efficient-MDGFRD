// author luigisbailo


double getR_aGF_nl ( struct particle_nl *particles, int *partList, double *distRow, int *distLabel,
                     double maxShell, int N, double L ){
//shell radius computed considering the expectd exit position of the other particles
//since it is considering only the exit positions, and not the actual shells it might overlap other shells
//this function is called after checking that the particle has space to construct a domain

  int iPart,jPart;
  double Dij[2];
  double Rshell[3],Rmin,R;
  double Shell;
  int j = 0;  

  iPart = particles[partList[0]].label;
  Rmin = particles[iPart].R_bd; 
  R = maxShell;

  //Dij is an array of two elements, the first-[0] contains the diff coeff of our particle,
    // the second-[1] contains the diff coeff of neighbour particles
  Dij[0] = particles[iPart].Diff;

  //it cycles over all particles determining the shell radius predicted according to the aGF algorithm
  //the lowest radius is finally selected
  while ( !(distLabel[j]<0) ) {


    jPart = distLabel[j];

    Dij[1] = particles[jPart].Diff;


    if ( particles[jPart].gf == false ) {


      if ( particles[jPart].type == particles[iPart].type )
        Shell = distRow [j] / 2;  //-particles[jPart].shellRed/2;
      else
        Shell = distRow [j] / (1 + particles[jPart].sqrtDiff/particles[iPart].sqrtDiff );

    }
    else{

      Rshell[0] = distRow[j] - particles[jPart].shell;

      Rshell[2] = sqrt(dist2next_per_nl(&particles[iPart],&particles[jPart],L))
                  - particles[iPart].radius - particles[jPart].radius;

      if ( particles[jPart].tau_exit - particles[iPart].time > Rshell[0]*Rshell[0]/6/Dij[0] ){  
      //if the first exit time of the j-part is higher than the average first exit time of the i-part from the largest domains,
      //the largest domain to i is proposed

        Rshell[1] = Rshell[0]; //- particles[iPart].shellRed;

      }
      else if ( Dij[0] != Dij[1] ) { 

        double deltaRoot = (Rshell[2]*Rshell[2])/(36*Dij[1]*Dij[1]) - 
                          ((1/(6*Dij[1])-1/(6*Dij[0]))*((Rshell[2]*Rshell[2]/
                                  (6*Dij[1]))+particles[jPart].tau_exit - particles[iPart].time));
        Rshell[1] = ( Rshell[2]/(6*Dij[1]) - sqrt(deltaRoot) ) / (1/(6*Dij[1])-1/(6*Dij[0]));
        // Rshell[1] -= particles[iPart].shellRed;
        
      }      
      else {

        Rshell[1] = Rshell[2]/2 + (3 * Dij[0] * (particles[jPart].tau_exit - particles[iPart].time) ) / Rshell[2];       
        // Rshell[1] -= particles[iPart].shellRed;

      } 

      Shell = fmin (Rshell[0],Rshell[1]);

    }
    Shell = Shell - particles[iPart].shellRed/2;
    if (Shell<R)  R=Shell;

    if (R<Rmin) {
      R=0;
      break;
    }
    
    j++;

  }

  return R;

}



double getR_GF_nl ( struct particle_nl *particles, int *partList, double *distRow, int *distLabel,
                    double maxShell, int N, double L ){

  int iPart,jPart;
  double Rshell;
  double minShell;
  double R = maxShell;
  int j = 0;

  iPart = particles[partList[0]].label;
  R = maxShell;

  if ( particles[iPart].gf == true )  
    minShell = particles[iPart].R_bd; 
  else 
    minShell = particles[iPart].R_gfrd;

  while ( !(distLabel[j]<0) ) {

    jPart = distLabel [j];
    // if (iPart==jPart) continue;

    // Rshell is an array of 3 elements, the 1st-[0] contains the distance with the j-particle,
      // the 2nd the distance with the expected exit position, the 3rd
    if ( particles[jPart].gf == false ){
      Rshell =  distRow [j] / 2;
    }
    else {
      Rshell = distRow [j] - particles[jPart].shell;
    }

    if( Rshell < minShell )
      return 0;
    else if (Rshell<R){
      R=Rshell;
    }

    j++;
  }

  return R;

}