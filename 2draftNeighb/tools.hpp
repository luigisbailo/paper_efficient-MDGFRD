class particle {
public:
  int label;
  int type;
  double pos [3];
  double pos_init [3];
  int pos_period [3];
  int idxBox [3];
  double shell;
  double time;
  double tau_exit;
  double pos_exit [3];
  double radius;
  double Diff;
  double sqrtDiff;
  double R_bd;
  double R_gfrd;
  double burstR;
  double shellRed;
  bool burst; 
  bool gf;
  double distBM;
  double distGF;
  double distBurst;
 };


  class BFdistances {
  public: 
    double dd[10];

  };


double dist2_per (particle *A, particle *B, double L ) {
//squared particles distance in periodic boundary conditions
 
 double XYZ [3];

  XYZ[0] = A->pos[0]-B->pos[0];
  if ( abs(XYZ[0]) > L/2 ) XYZ[0] = L-abs(XYZ[0]);
  XYZ[1] = A->pos[1]-B->pos[1];
  if ( abs(XYZ[1]) > L/2 ) XYZ[1] = L-abs(XYZ[1]);
  XYZ[2] = A->pos[2]-B->pos[2];
  if ( abs(XYZ[2]) > L/2 ) XYZ[2] = L-abs(XYZ[2]);

  return XYZ[0]*XYZ[0] + XYZ[1]*XYZ[1] + XYZ[2]*XYZ[2];  

}

double dist2next_per (particle *A, particle *B, double L ) {

  double XYZ [3];

  XYZ[0] = A->pos[0]-B->pos_exit[0];
  if ( abs(XYZ[0]) > L/2 ) XYZ[0] = L-abs(XYZ[0]);
  XYZ[1] = A->pos[1]-B->pos_exit[1];
  if ( abs(XYZ[1]) > L/2 ) XYZ[1] = L-abs(XYZ[1]);
  XYZ[2] = A->pos[2]-B->pos_exit[2];
  if ( abs(XYZ[2]) > L/2 ) XYZ[2] = L-abs(XYZ[2]);

  return XYZ[0]*XYZ[0] + XYZ[1]*XYZ[1] + XYZ[2]*XYZ[2];

}


void getDist ( particle *particles, vector<vector<int>> *grid, double *distRow, int *distLabel, int ilabel, int N, double L, int nBoxes ) {
  //The j-th element of the vector "distRow[]" is filled with the distance with the j-th particle
  //The j-th particle means the j-th ordered on a exit-time basis
  
  int pos;
  int countDist = 0;
  int iBoxi =  particles[ilabel].idxBox[0];
  int iBoxj =  particles[ilabel].idxBox[1];
  int iBoxk =  particles[ilabel].idxBox[2];
  int i,j,k;

  while (!(distLabel[countDist]<0)){
    distRow [countDist] = -1;
    distLabel [countDist] = -1;
    countDist ++;
  }

// cout << nBoxes<<endl;


  countDist=0;

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

        // cout << i << "\t" << j << "\t" << k << endl;
        for(int nPart=1; nPart<(*grid)[i*nBoxes*nBoxes+j*nBoxes+k][0]+1; nPart++ ){
          
          int jlabel = (*grid)[i*nBoxes*nBoxes+j*nBoxes+k][nPart];

          if (jlabel==ilabel) continue;

          distLabel[countDist] = jlabel;

          if (particles[jlabel].gf==true){ 
            distRow [countDist] = sqrt(dist2_per ( &particles[ilabel], &particles[jlabel], L )) - particles[ilabel].radius - particles[jlabel].radius;
          }   
          else if ( abs(particles[ilabel].time-particles[jlabel].time) < 0.00001 ){
            distRow [countDist] = sqrt(dist2_per ( &particles[ilabel], &particles[jlabel], L )) - particles[ilabel].radius - particles[jlabel].radius;
          }  
          else if ( abs(particles[ilabel].time-particles[jlabel].tau_exit) < 0.00001 ){
            distRow [countDist] = sqrt(dist2next_per ( &particles[ilabel], &particles[jlabel], L )) - particles[ilabel].radius - particles[jlabel].radius;
          }
          else {
            cout << "ERROR: distance in getDist" << endl;
          }

          countDist++;

        }
      }
    }
  }


}



bool compareTime (particle A, particle B)  {return A.tau_exit < B.tau_exit;}


void sortPart ( particle *particles, int *partList, int N) {

  int tempList;

  for (int n=0; n<N-1; n++){

    if (particles[partList[n]].tau_exit>particles[partList[n+1]].tau_exit){

      tempList=partList[n];
      partList[n]=partList[n+1];
      partList[n+1]=tempList;

    }
    else break;

  }

}

void printPos_per (particle *particle, int *partList, int N);

void sortBurst ( int burstLabel, int *partList, int N) {

  int tempList,n;

  for ( n=1; n<N; n++){

    if ( burstLabel== partList[n] )
      break;

  }

  tempList = partList[n];

  for ( int m = n; m>1; m-- ){

    partList[m] = partList[m-1];

  }
  partList[1]=tempList;


}



void checkBound (double *pos, int *pos_period, double L) {

  if ( pos[0] > L ) {
    pos[0] -= L;
    pos_period[0] ++;
  }
  else if ( pos[0] < 0 ) {
    pos[0] += L;
    pos_period[0] --;
  }

  if ( pos[1] > L ) {
    pos[1] -= L;
    pos_period[1] ++;
  }
  else if ( pos[1] < 0 ) {
    pos[1] += L;
    pos_period[1] --;
  }

  if ( pos[2] > L ) {
    pos[2] -= L;
    pos_period[2] ++;
  }
  else if ( pos[2] < 0 ) {
    pos[2] += L;
    pos_period[2] --;
  }

};


void fixBound_aGF ( double *pos, double *pos_exit, int *pos_period, double L) {

  if ( pos[0]-pos_exit[0] > L/2 )
    pos_period [0] --;
  else  if ( pos_exit[0]-pos[0] > L/2 )
      pos_period[0] ++;

  if ( pos[1]-pos_exit[1] > L/2 )
    pos_period [1] --;
  else if ( pos_exit[1]-pos[1] > L/2 )
      pos_period[1] ++;

  if ( pos[2]-pos_exit[2] > L/2 )
    pos_period [2] --;
  else if ( pos_exit[2]-pos[2] > L/2 )
      pos_period[2] ++;

};


void polarTransf ( double *pos, double R, double u, double v ) {
//Theta is defined in [0,2pi]
//Phi is defined in [0,pi]
  
  double theta = 2 * M_PI * u;
  double phi = acos( 2*v - 1 );

  pos[0] = R * cos(theta) * sin(phi);      
  pos[1] = R * sin(theta) * sin(phi);
  pos[2] = R * cos(phi);  

}


void updatePart_aGF ( particle *P, vector <vector <int>> *grid, gsl_rng *r, double dt, double L, double boxSizeGF, int nBoxesGF ) {
// updates are already performed in case of bursting

  if ( P->burst == false) {


    int i,j,k,iex,jex,kex;
    i = P->idxBox[0];
    j = P->idxBox[1];
    k = P->idxBox[2];
    iex = trunc ( P->pos_exit[0]/boxSizeGF );
    jex = trunc ( P->pos_exit[1]/boxSizeGF );
    kex = trunc ( P->pos_exit[2]/boxSizeGF );
    if ( ( (i==iex)&&(j==jex)&&(k==kex) ) == false ){
      
      int label = P -> label;

      //remove particle from previous grid
      int nElem = (*grid)[i*nBoxesGF*nBoxesGF + j*nBoxesGF +k][0]; 
      int n;
      for ( n=1; n<nElem+1; n++ ){
        if ( label == (*grid)[i*nBoxesGF*nBoxesGF + j*nBoxesGF +k][n]){
          (*grid)[i*nBoxesGF*nBoxesGF + j*nBoxesGF +k][n] = -1;
          break;
        }
      }
      while (n<nElem){
        (*grid)[i*nBoxesGF*nBoxesGF + j*nBoxesGF +k][n]=(*grid)[i*nBoxesGF*nBoxesGF + j*nBoxesGF +k][n+1];
        (*grid)[i*nBoxesGF*nBoxesGF + j*nBoxesGF +k][n+1] = -1;
        n++;
      }
      (*grid)[i*nBoxesGF*nBoxesGF + j*nBoxesGF +k][0] --;

    
      //add particle to new grid
      (*grid)[iex*nBoxesGF*nBoxesGF + jex*nBoxesGF +kex][0]++;
      nElem = (*grid)[iex*nBoxesGF*nBoxesGF + jex*nBoxesGF +kex][0];
      (*grid)[iex*nBoxesGF*nBoxesGF + jex*nBoxesGF +kex][nElem]=label;


      //update particle's indexes
      P ->idxBox[0] = iex;
      P ->idxBox[1] = jex;
      P ->idxBox[2] = kex;


    }

    P->pos[0] = P->pos_exit[0];
    P->pos[1] = P->pos_exit[1];
    P->pos[2] = P->pos_exit[2];
    P->shell = 0;
    P->time = P->tau_exit;

  }
  if ( P->gf == true){
  
    double deltaT = trunc( (P->time+dt)/dt )*dt - P->time;

    P -> pos[0] += gsl_ran_gaussian (r,1)*P->sqrtDiff*sqrt(2*deltaT);
    P -> pos[1] += gsl_ran_gaussian (r,1)*P->sqrtDiff*sqrt(2*deltaT);
    P -> pos[2] += gsl_ran_gaussian (r,1)*P->sqrtDiff*sqrt(2*deltaT);
    checkBound ( P->pos, P->pos_period, L );
    P->pos_exit[0] = P->pos[0];
    P->pos_exit[1] = P->pos[1];
    P->pos_exit[2] = P->pos[2];
    P->time += deltaT;
    P->tau_exit = P->time;

  }


}



void updatePart_GF ( particle *P, vector<vector<int>> *grid, gsl_rng *r, double dt, double L, double boxSizeGF, int nBoxesGF ) {
// differently from aGF positions are not sampled when the domain is constructed

  double deltaPos [3];

  if ( P->gf == true ){

    polarTransf ( deltaPos, P -> shell, gsl_rng_uniform (r), gsl_rng_uniform (r));
    //deltaPos now contains the displacements in cartesian coordinates
    P ->pos_exit[0] = P -> pos[0] + deltaPos[0];  
    P ->pos_exit[1] = P -> pos[1] + deltaPos[1];
    P ->pos_exit[2] = P -> pos[2] + deltaPos[2];      
    P->shell = 0;
    P->time = P->tau_exit;

    double deltaT = trunc( (P->time+dt)/dt )*dt - P->time;

    P -> pos_exit[0] += gsl_ran_gaussian (r,1)*P->sqrtDiff*sqrt(2*deltaT);
    P -> pos_exit[1] += gsl_ran_gaussian (r,1)*P->sqrtDiff*sqrt(2*deltaT);
    P -> pos_exit[2] += gsl_ran_gaussian (r,1)*P->sqrtDiff*sqrt(2*deltaT);
    checkBound ( P->pos_exit, P->pos_period, L );
    P->time += deltaT;
    P->tau_exit = P->time;

    int i,j,k,iex,jex,kex;
    i = P->idxBox[0];
    j = P->idxBox[1];
    k = P->idxBox[2];
    iex = trunc ( P->pos_exit[0]/boxSizeGF );
    jex = trunc ( P->pos_exit[1]/boxSizeGF );
    kex = trunc ( P->pos_exit[2]/boxSizeGF );
    if ( ( (i==iex)&&(j==jex)&&(k==kex) ) == false ){
      
      int label = P -> label;

      //remove particle from previous grid
      int nElem = (*grid)[i*nBoxesGF*nBoxesGF + j*nBoxesGF +k][0]; 
      int n;
      for ( n=1; n<nElem+1; n++ ){
        if ( label == (*grid)[i*nBoxesGF*nBoxesGF + j*nBoxesGF +k][n]){
          (*grid)[i*nBoxesGF*nBoxesGF + j*nBoxesGF +k][n] = -1;
          break;
        }
      }
      while (n<nElem){
        (*grid)[i*nBoxesGF*nBoxesGF + j*nBoxesGF +k][n]=(*grid)[i*nBoxesGF*nBoxesGF + j*nBoxesGF +k][n+1];
        (*grid)[i*nBoxesGF*nBoxesGF + j*nBoxesGF +k][n+1] = -1;
        n++;
      }
      (*grid)[i*nBoxesGF*nBoxesGF + j*nBoxesGF +k][0] --;

    
      //add particle to new grid
      (*grid)[iex*nBoxesGF*nBoxesGF + jex*nBoxesGF +kex][0]++;
      nElem = (*grid)[iex*nBoxesGF*nBoxesGF + jex*nBoxesGF +kex][0];
      (*grid)[iex*nBoxesGF*nBoxesGF + jex*nBoxesGF +kex][nElem]=label;


      //update particle's indexes
      P ->idxBox[0] = iex;
      P ->idxBox[1] = jex;
      P ->idxBox[2] = kex;

    }

    P->pos[0] = P->pos_exit[0];
    P->pos[1] = P->pos_exit[1];
    P->pos[2] = P->pos_exit[2];

  }
  else if ( P->burst == false ){
    // it is the case of a BM integration

    int i,j,k,iex,jex,kex;
    i = P->idxBox[0];
    j = P->idxBox[1];
    k = P->idxBox[2];
    iex = trunc ( P->pos_exit[0]/boxSizeGF );
    jex = trunc ( P->pos_exit[1]/boxSizeGF );
    kex = trunc ( P->pos_exit[2]/boxSizeGF );
    if ( ( (i==iex)&&(j==jex)&&(k==kex) ) == false ){
      
      int label = P -> label;

      //remove particle from previous grid
      int nElem = (*grid)[i*nBoxesGF*nBoxesGF + j*nBoxesGF +k][0]; 
      int n;
      for ( n=1; n<nElem+1; n++ ){
        if ( label == (*grid)[i*nBoxesGF*nBoxesGF + j*nBoxesGF +k][n]){
          (*grid)[i*nBoxesGF*nBoxesGF + j*nBoxesGF +k][n] = -1;
          break;
        }
      }
      while (n<nElem){
        (*grid)[i*nBoxesGF*nBoxesGF + j*nBoxesGF +k][n]=(*grid)[i*nBoxesGF*nBoxesGF + j*nBoxesGF +k][n+1];
        (*grid)[i*nBoxesGF*nBoxesGF + j*nBoxesGF +k][n+1] = -1;
        n++;
      }
      (*grid)[i*nBoxesGF*nBoxesGF + j*nBoxesGF +k][0] --;

    
      //add particle to new grid
      (*grid)[iex*nBoxesGF*nBoxesGF + jex*nBoxesGF +kex][0]++;
      nElem = (*grid)[iex*nBoxesGF*nBoxesGF + jex*nBoxesGF +kex][0];
      (*grid)[iex*nBoxesGF*nBoxesGF + jex*nBoxesGF +kex][nElem]=label;


      //update particle's indexes
      P ->idxBox[0] = iex;
      P ->idxBox[1] = jex;
      P ->idxBox[2] = kex;


    }

    P->pos[0] = P->pos_exit[0];
    P->pos[1] = P->pos_exit[1];
    P->pos[2] = P->pos_exit[2];
    P->time = P->tau_exit;

  }

}


void relabel_part (particle *particles, int N){

  for ( int n=0; n<N; n++){
    particles[n].label = n;
  }


}