// double drawTimeNewt ( double b, double D, double xi ) {

//   double t,tmem;

//   t = 0.0917517*b*b/D;
//   tmem = -1;

//   double S = 1-Sfunct (t,b,D);
//   double dS = Sder (t,b,D);
//   int count = 0;


//   while ( abs(t-tmem) > 0.0001 | abs(S-xi) > DRAW_CONVERGENCE ) {

    

//     // cout << count << "\t" << t << "\t" << S-xi << "\t" << (S-xi)/dS << endl;

//     count++;
//     if (count > MAX_ITERATIONS){
//       cout << setprecision (15);
//       cout << "Error: finding the root of S was not possible for:" << endl;
//       cout << "b = " << b << "\t D = " << D << "\t xi = " << xi << endl;     
//       return t;
//     }

//     tmem = t;
//     t = t + (S-xi)/dS;

//     S = 1-Sfunct (t,b,D);
//     if ( S==1 ) return t;
//     dS = Sder (t,b,D);

//   }

//   return t;
    
// }


// double drawPosNewt ( double t, double b, double D, double xi ) {

//   double r;

//   double t0=0.063;
//   double t1=0.234;
//   if (t<t0*b*b/D) r = sqrt(t*D)*2;
//   else if (t<t1*b*b/D) {
//     double R0=2*sqrt(t0*b*b);
//     double R1=0.646*b;
//     double beta = (R0*exp(pow(t0,0.5)+pow(t1,0.5))-R1*exp(pow(t0,0.5)+pow(t1,0.5)))/(exp(pow(t0,0.5))-exp(pow(t1,0.5)));
//     double gamma = -(R0*(exp(pow(t1,0.5))-1)*exp(pow(t0,0.5))-R1*(exp(pow(t0,0.5))-1)*exp(pow(t1,0.5)))  / (exp(pow(t0,0.5))-exp(pow(t1,0.5)));
//     r=beta*(1-exp(-pow(t*D/b/b,0.5)))+gamma;
//   }
//   else r = 0.646*b;


//     double S = Sfunct (t,b,D);
//     double P = Pfunct (r,t,b,D,S);
//     double dP = Pder (r,t,b,D,S);
//     int count = 0;
//     double rMem=-1;

 
//     while ( abs(r-rMem) > 0.0001 | abs(P-xi) > DRAW_CONVERGENCE ) {

//       count++;
//       if (count > MAX_ITERATIONS){
//         cout << setprecision (15);
//         cout << "Error: finding the root of P was not possible for:" << endl;
//         cout << "t = " << t << "\tb = " << b << "\t D = " << D << "\t xi = " << xi << endl;     
//         return r;
//       }

//       rMem = r;
//       r = r - (P-xi)/dP;

//       P = Pfunct (r,t,b,D,S);
//       dP = Pder (r,t,b,D,S);

//     }

//   return r;

// }

