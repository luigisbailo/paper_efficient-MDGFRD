//TO COMPILE: g++ -std=c++11 main.cpp -o main -lgsl -lgslcblas -lm

using namespace std;
#define print(x) cout << x << endl;

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_math.h>
#include <cmath>
#include <gsl/gsl_errno.h>
#include <iostream>
#include <fstream>
#include <algorithm> 
#include <vector>
#include <iomanip>
#include <chrono>
#include <cstring>
#include <sstream>

using namespace std::chrono;

#include "parameters.hpp"
#include "tools.hpp"
#include "greensFunct.hpp"
#include "draw.hpp"
#include "init.hpp"
#include "step.hpp"
#include "shell.hpp"
#include "print.hpp"
#include "burst.hpp"
#include "bruteForce.hpp"
#include "checks.hpp"

int main () {         

cout << drawTimeNewt (233.96638131678,0.01,0.999997306382284) << endl;
// cout << drawPosNewt (233.96638131678,80.0079415409098,0.01,0.992160857422277) << endl;

  }




