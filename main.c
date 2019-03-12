#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_math.h>

#include "LowParticleNumber/parameters.h"
#include "LowParticleNumber/tools.h"
#include "LowParticleNumber/bruteForce.h"
#include "LowParticleNumber/greensFunct.h"
#include "LowParticleNumber/draw.h"
#include "LowParticleNumber/init.h"
#include "LowParticleNumber/burst.h"
#include "LowParticleNumber/shell.h"
#include "LowParticleNumber/step.h"
#include "LowParticleNumber/print.h"
#include "LowParticleNumber/checks.h"
#include "LowParticleNumber/run_aGF1.h"
#include "LowParticleNumber/run_aGF2.h"
#include "LowParticleNumber/run_BM.h"
#include "LowParticleNumber/run_GF.h"
#include "LowParticleNumber/run_hybrGF.h"
#include "LowParticleNumber/Fig4/Fig4.h"
#include "LowParticleNumber/Fig5_8/Fig5_8.h"
#include "LowParticleNumber/Fig9/Fig9.h"
#include "LowParticleNumber/alpha/fixAlpha/fixAlpha.h"

#include "NeighborLists/parameters_nl.h"
#include "NeighborLists/tools_nl.h"
#include "NeighborLists/bruteForce_nl.h"
#include "NeighborLists/init_nl.h"
#include "NeighborLists/burst_nl.hpp"
#include "NeighborLists/shell_nl.h"
#include "NeighborLists/step_nl.h"
#include "NeighborLists/print_nl.h"
#include "NeighborLists/checks_nl.h"
#include "NeighborLists/run_GF_nl.h"
#include "NeighborLists/run_aGF_nl.h"
#include "NeighborLists/run_BM_nl.h"
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

    printf("\n\nFig9b is being produced:\n\n");

    D_A = 0.01;
    D_B = 0.001;
    R_A = 1.5;
    R_B = 3.5;

    fig9 ( D_A, D_B, R_A, R_B );

    printf("\n\nTab1 is being produced:\n\n");

    tab1(210);

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