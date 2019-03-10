
//#include <math.h>
//#include <stdio.h>
//#include <stdlib.h>
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
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_math.h>

#include "LowParticleNumber/parameters.hpp"
#include "LowParticleNumber/tools.hpp"
#include "LowParticleNumber/bruteForce.hpp"
#include "LowParticleNumber/greensFunct.hpp"
#include "LowParticleNumber/draw.hpp"
#include "LowParticleNumber/init.hpp"
#include "LowParticleNumber/burst.hpp"
#include "LowParticleNumber/shell.hpp"
#include "LowParticleNumber/step.hpp"
#include "LowParticleNumber/print.hpp"
#include "LowParticleNumber/checks.hpp"
#include "LowParticleNumber/run_aGF1.hpp"
#include "LowParticleNumber/run_aGF2.hpp"
#include "LowParticleNumber/run_BM.hpp"
#include "LowParticleNumber/run_GF.hpp"
#include "LowParticleNumber/run_hybrGF.hpp"
#include "LowParticleNumber/Fig4/Fig4.hpp"
#include "LowParticleNumber/Fig5_8/Fig5_8.hpp"
#include "LowParticleNumber/Fig9/Fig9.hpp"
#include "LowParticleNumber/alpha/fixAlpha/fixAlpha.hpp"

#include "NeighborLists/parameters.hpp"
#include "NeighborLists/tools_nl.hpp"
#include "NeighborLists/bruteForce_nl.hpp"
#include "NeighborLists/init_nl.hpp"
#include "NeighborLists/burst_nl.hpp"
#include "NeighborLists/shell_nl.hpp"
#include "NeighborLists/step_nl.hpp"
#include "NeighborLists/print.hpp"
#include "NeighborLists/checks_nl.hpp"
#include "NeighborLists/run_GF_nl.hpp"
#include "NeighborLists/run_aGF_nl.hpp"
#include "NeighborLists/run_BM_nl.hpp"
#include "NeighborLists/Tab1/tab1.hpp"



int main () {


    double D_A, D_B, R_A,R_B;

    printf("\n\nFig4 is being produced:\n\n");

    fig4 ();

    printf("\n\nFig5-8a are being produced:\n\n");

    D_A = 0.01;
    D_B = 0.01;
    R_A = 2.5;
    R_B = 2.5;

    fig5_8 ( D_A, D_B, R_A, R_B );

    printf("\n\nFig5-8b are being produced:\n\n");

    D_A = 0.01;
    D_B = 0.001;
    R_A = 1.5;
    R_B = 3.5;

    fig5_8 ( D_A, D_B,  R_A, R_B  );

    printf("\n\nFig9a is being produced:\n\n");

    D_A = 0.01;
    D_B = 0.01;
    R_A = 2.5;
    R_B = 2.5;

    fig9 ( D_A, D_B, R_A, R_B );

    D_A = 0.01;
    D_B = 0.001;
    R_A = 1.5;
    R_B = 3.5;

    printf("\n\nFig9b is being produced:\n\n");

    fig9 ( D_A, D_B, R_A, R_B );

    printf("\n\nTab1 is being produced:\n\n");

    tab1(210);
    tab1(2100);

    printf("\n\nAppendix is being produced:\n\n");

    double L;

    L =  1150;
    fixAlpha(L);

    L = 545;
    fixAlpha(L);

    L = 255;
    fixAlpha(L);

    L = 118;
    fixAlpha(L);

    L = 55;
    fixAlpha(L);

    L = 25;
    fixAlpha(L);



}