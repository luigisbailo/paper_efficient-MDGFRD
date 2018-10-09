void printPos ( particle *particles, int *partList, int N){

  cout <<"Part.\t"<<"Time\t"<<"Tau\t"<<"Shell\t"
       <<"x\t"<<"y\t"<<"z\t"
       <<"x_ex\t"<<"y_ex\t"<<"z_ex\t"
       <<"Radius\t" <<"BURST\n";

  for (int count=0; count<N; count++){
    cout <<particles[partList[count]].label << "\t"
	 <<particles[partList[count]].time << "\t" 
	 <<particles[partList[count]].tau_exit<<"\t"
	 <<particles[partList[count]].shell<<"\t"
	 <<particles[partList[count]].pos[0]<<"\t"
	 <<particles[partList[count]].pos[1]<<"\t"
	 <<particles[partList[count]].pos[2] <<"\t"
	 <<particles[partList[count]].pos_exit[0]<<"\t"
	 <<particles[partList[count]].pos_exit[1]<<"\t"
	 <<particles[partList[count]].pos_exit[2] <<"\t"
	 <<particles[partList[count]].radius <<"\t"
	 <<particles[partList[count]].burst <<"\n";
  }
}

// void printPos_per (particle *particles, int *partList, int N){

//   cout <<"Part.\t"<<"Time\t\t"<<"Tau\t\t"<<"Shell\t"
//        <<"x\t"<<"y\t"<<"z\t"
//        <<"x_ex\t"<<"y_ex\t"<<"z_ex\t"
//        <<"x_per\t"<<"y_per\t"<<"z_per\t"
//        <<"Radius\t" <<"BURST\t"<<"GF\n";

//   for (int count=0; count<N; count++){
//     cout <<particles[partList[count]].label << "\t"
// 	 <<particles[partList[count]].time << "\t\t" 
// 	 <<particles[partList[count]].tau_exit<<"\t\t"
// 	 <<particles[partList[count]].shell<<"\t"
// 	 <<particles[partList[count]].pos[0]<<"\t"
// 	 <<particles[partList[count]].pos[1]<<"\t"
// 	 <<particles[partList[count]].pos[2] <<"\t"
// 	 <<particles[partList[count]].pos_exit[0]<<"\t"
// 	 <<particles[partList[count]].pos_exit[1]<<"\t"
// 	 <<particles[partList[count]].pos_exit[2] <<"\t"
// 	 <<particles[partList[count]].pos_period[0]<<"\t"
// 	 <<particles[partList[count]].pos_period[1]<<"\t"
// 	 <<particles[partList[count]].pos_period[2] <<"\t"
// 	 <<particles[partList[count]].radius <<"\t"
// 	 <<particles[partList[count]].burst <<"\t"
// 	 <<particles[partList[count]].gf<<"\n";
//   }
// }


void printDist_per (particle *particles, int *partList, int N, double L ) {
 
  for (int n=1; n<N; n++){
	double dist = sqrt(dist2_per(&particles[partList[0]],&particles[partList[n]],L));
	double distNext = sqrt(dist2next_per(&particles[partList[0]],&particles[partList[n]],L));
	double ishell = particles[partList[0]].shell;
	double jshell = particles[partList[n]].shell;
	cout << "R "<<particles[partList[0]].label<< "-" <<particles[partList[n]].label<<" = "<< dist << "\t"
	     << "R - shells -inter = "<<dist - (ishell+jshell+particles[partList[n]].radius+particles[partList[0]].radius) << "\t"
	     << "Rnext -inter = "<<distNext - (particles[partList[n]].radius+particles[partList[0]].radius) 
	     << "\n";  

  }

}



void print_shell ( particle *particles, int N){
  
  for (int i=0; i<N; i++){
    cout << "Particle " << i << "\tDiff = " << particles[i].Diff 
	 << "\tR_bm = " << particles[i].R_bd
	 << "\tShell = " << particles[i].shell 
	 << "\tTau = " << particles[i].tau_exit << "\n";
            
  }
}


