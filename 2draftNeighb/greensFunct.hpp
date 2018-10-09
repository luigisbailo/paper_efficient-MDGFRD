// double Sfunct ( double t, double b, double D) {
// //This function returns 1-Survival probability, i.e. the probability of having already left the domain

//   double S = 0;
//   double term,termA,termB;
//   int m = 1;
//   double conv = 0.000000000000000000001;

//   if ( t>= b*b/D/1000 ) {

//     do {

//       termA = exp ( -m*m*M_PI*M_PI*D*t/(b*b) );
//       termB = exp ( -(m+1)*(m+1)*M_PI*M_PI*D*t/(b*b) );
//       term = 2 * (termA-termB);


//       S += term;
//       m += 2;

//     } while ( abs (term) > conv | S < 0 );

//     S = 1-S;
//   }
//   else 
//     S = 0;

//   return S;

// }



// double Sder ( double t, double b, double D) {

//   double S=0;
//   double term,termA,termB;
//   int m=1;
//   double conv = 0.0000000000000000000001/(b*b/D);


//   if ( t>= b*b/D/1000 ) {

//     do {

//       termA = m * m * exp ( - m*m*M_PI*M_PI*D*t/(b*b) );
//       termB = (m+1) * (m+1) * exp ( - (m+1)*(m+1)*M_PI*M_PI*D*t/(b*b) );
//       term = 2*D*M_PI*M_PI/b/b * ( termA - termB );

//       S += term;
//       m += 2;

//     } while ( abs (term) > conv | S < 0  );
//   }
//   else 
//     S = 0;

//   return S;

// }


// double Pfunct ( double radius, double t, double b, double D, double S ) {

//   double P=0;
//   double  term, termA, termB;
//   int m=1;
//   double conv = 0.00000000000000001;

      
//   do {

//     termA = exp ( - M_PI*M_PI*D*t*m*m/(b*b));
//     termB =  b*sin(m*M_PI*radius/b)/m - M_PI*radius*cos(m*M_PI*radius/b); 
//     term =  2/(b*M_PI) * termA * termB / (1-S);
//     P += term;
//     m += 1;

//   } while ( abs(term) > conv  |  termA>termB/10000000);


//   return P;  

// }


// double Pder ( double radius, double t, double b, double D, double S ) {

//   double P=0;
//   double term,termA,termB;
//   int m=1;
//   double conv = 0.00000000000000001*b;


//   do {

//     termA = exp ( -M_PI*M_PI*D*t*m*m/(b*b)  );
//     termB = sin (m*M_PI*radius/b) * m * radius;
//     term = 2*M_PI/b/b * termA * termB / (1-S);
//     P += term;
//     m += 1;

//   } while ( abs(term) > conv  | termA>termB/10000000 );

//   return P;
// }
