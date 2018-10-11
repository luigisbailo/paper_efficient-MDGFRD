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
#include "run_aGF1.hpp"
#include "run_aGF2.hpp"
#include "run_BM.hpp"
#include "run_GF.hpp"
#include "run_hybrGF.hpp"
#include "alpha/GF-BF/GF_BF.hpp"
#include "performance/performance.hpp"
#include "diffusion/diffusion.hpp"
#include "alpha/fixAlpha/fixAlpha.hpp"



int main () {         


    std::cout << "Fig. 4" << std::endl<<std::endl;
    GF_BF();
    std::cout << std::endl << std::endl;


    std::cout << "Fig. 5, Fig. 6, Fig.7, Fig. 8" << std::endl<<std::endl;
    performance();
    std::cout << std::endl << std::endl;


    std::cout << "Fig. 9" << std::endl<<std::endl;
    diffusion();
    std::cout << std::endl << std::endl;


    std::cout << "Fig. 10" << std::endl<<std::endl;
    fixAlpha();
    std::cout << std::endl << std::endl;


  }




