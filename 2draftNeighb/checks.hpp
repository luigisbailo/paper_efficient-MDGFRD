void check_aGF (  particle *particles, int *partList, int N, double L ){
  

  if (particles[partList[0]].time > particles[partList[0]].tau_exit){
    cout << "ERROR: exit time\n";
    // exit (EXIT_FAILURE);
  }

  if ( particles[partList[0]].pos[0]>L | particles[partList[0]].pos[0]<0 | particles[partList[0]].pos[1]>L | particles[partList[0]].pos[1] <0  | particles[partList[0]].pos[2]>L | particles[partList[0]].pos[2]<0  ) {
    cout << "ERROR: pos\n";
    // printPos (particles,N);
    // exit (EXIT_FAILURE);
  
  } 

  if (  particles[partList[0]].gf == true  &&
       (particles[partList[0]].pos_exit[0]>L | particles[partList[0]].pos_exit[0]<0 | particles[partList[0]].pos_exit[1]>L | particles[partList[0]].pos_exit[1] <0  | particles[partList[0]].pos_exit[2]>L | particles[partList[0]].pos_exit[2]<0 ) ) {
    cout << "ERROR: pos_exit\n";
    // printPos (particles,N);
    // exit (EXIT_FAILURE);
  
  } 

  for (int count=1; count<N; count++){
    int jPart = partList[count];
    if ( sqrt(dist2_per(&particles[partList[0]],&particles[jPart],L)) - particles[partList[0]].shell - particles[jPart].shell - particles[partList[0]].radius - particles[jPart].radius < - 0.0000000001  && particles[jPart].shell>0 ){
      cout << "ERROR: distance particles " << particles[partList[0]].label << "-" <<  particles[jPart].label << "\n";
      cout << sqrt(dist2_per(&particles[partList[0]],&particles[jPart],L)) - particles[partList[0]].shell - particles[jPart].shell << "\n";
      // exit(EXIT_FAILURE);
    }
    
    // if (sqrt(dist2_per(&particles[partList[0]],&particles[jPart],XYZ))<0.001) {
     //   cout << "ERROR: repulsion\n";
     //   //exit(EXIT_FAILURE);
    // }
  }
  
}


void check_GF (  particle *particles, int *partList, int N, double L ){
  

  if (particles[partList[0]].time > particles[partList[0]].tau_exit){
    cout << "ERROR: exit time\n";
    exit (EXIT_FAILURE);
  }


  if ( particles[partList[0]].pos[0]>L | particles[partList[0]].pos[0]<0 | particles[partList[0]].pos[1]>L | particles[partList[0]].pos[1] <0  | particles[partList[0]].pos[2]>L | particles[partList[0]].pos[2]<0  ) {
    cout << "ERROR: pos\n";
    // printPos (particles,N);
    exit (EXIT_FAILURE);
  
  } 

  for (int count=1; count<N; count++){

    int jPart = partList[count];
    if ( sqrt(dist2_per(&particles[partList[0]],&particles[jPart],L)) - particles[partList[0]].shell - particles[jPart].shell - particles[partList[0]].radius - particles[jPart].radius < - 0.0000000001  && particles[jPart].shell>0 ){
      cout << "ERROR: distance particles " << particles[partList[0]].label << "-" <<  particles[jPart].label << "\n";
      cout << sqrt(dist2_per(&particles[partList[0]],&particles[jPart],L)) - particles[partList[0]].shell - particles[jPart].shell << "\n";
      exit(EXIT_FAILURE);
    }
    
    // if (sqrt(dist2_per(&particles[partList[0]],&particles[jPart],XYZ))<0.001) {
	   //   cout << "ERROR: repulsion\n";
	   //   //exit(EXIT_FAILURE);
    // }
  }
  
}


void check_times ( particle *particles, int *partList, int N) {

  for (int count=1; count<N; count++){

    if ( particles[count].tau_exit<particles[partList[0]].time){
      cout << "ERROR: time1\n";
      cout << particles[partList[0]].time << "\t" << particles[count].tau_exit << "\n";
      exit (EXIT_FAILURE);
    }
    if ( particles[partList[0]].time<particles[count].time){
      cout << "ERROR: time2\n";
      cout << particles[partList[0]].time << "\t" << particles[count].time << "\n";
    printPos (particles,partList,N);
      exit (EXIT_FAILURE);
    }

  }

}

void check_Neighb ( particle *particles, vector<vector<int>> *grid, int N, int nBoxes, double boxSize, int NlimBox ){

  int i,j,k;
  for ( int n=0; n<N; n++){
    bool present=false;
    i = particles[n].idxBox[0];
    j = particles[n].idxBox[1];
    k = particles[n].idxBox[2];
    int nPartBox = (*grid)[i*nBoxes*nBoxes + j*nBoxes +k][0];
    for ( int count=1; count<nPartBox+1; count++ ){
      if ( (*grid)[i*nBoxes*nBoxes + j*nBoxes +k][count] == n ) present = true;
    }
    if ( present==false ){
      cout << "ERROR: not present NeighbList" << endl;
      // exit(EXIT_FAILURE);
    }
  } 

  for ( i=0; i<nBoxes; i++ ){
    for ( j=0; j<nBoxes; j++ ){
      for ( k=0; k<nBoxes; k++ ){

        int nPartBox = (*grid)[i*nBoxes*nBoxes + j*nBoxes +k][0];
        int count = 0;
        while (count<nPartBox){
          count++;
          int label =  (*grid)[i*nBoxes*nBoxes + j*nBoxes +k][count];
          if ( !(i==particles[label].idxBox[0] && j==particles[label].idxBox[1] && k==particles[label].idxBox[2] ) ){
            cout << "ERROR: grid->particle NeighbList" << endl;
            exit(EXIT_FAILURE);
          }
        } 
        count++;
        while (count<NlimBox){
          if ( (*grid)[i*nBoxes*nBoxes + j*nBoxes +k][count] != -1 ) {
            cout << "ERROR: grid filling NeighbList" << endl;
            exit(EXIT_FAILURE);
          }
          count++;
        }

    
      }
    }
  }


}


// void checkShells (vector <particle> *particles, int N, double *XYZ){

//   if (particles[0].shell>0) {
//     for (int j=1; j<N; j++) {

//       if (sqrt(dist2_per(&particles[0],&particles[j],XYZ)) < particles[0].shell + particles[j].shell + particles[0].radius + particles[j].radius) {
//         cout << "ERROR: shells \t" 
//              << sqrt(dist2_per(&particles[0],&particles[j],XYZ)) - particles[0].shell - particles[j].shell - particles[0].radius - particles[j].radius
//              << "\n"; 
//         // exit(EXIT_FAILURE);
//       }
//     }

//   }


// }
