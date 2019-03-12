cmake_minimum_required(VERSION 3.9)
project(efficient_MDGFRD)

set(CMAKE_C_STANDARD 11)

add_executable(efficient_MDGFRD
        main.c
        LowParticleNumber/Fig4/Fig4.h
        LowParticleNumber/Fig9/Fig9.h
        LowParticleNumber/bruteForce.h
        LowParticleNumber/burst.h
        LowParticleNumber/checks.h
        LowParticleNumber/draw.h
        LowParticleNumber/greensFunct.h
        LowParticleNumber/init.h
        LowParticleNumber/parameters.h
        LowParticleNumber/print.h
        LowParticleNumber/run_aGF1.h
        LowParticleNumber/run_aGF2.h
        LowParticleNumber/run_BM.h
        LowParticleNumber/run_GF.h
        LowParticleNumber/run_hybrGF.h
        LowParticleNumber/shell.h
        LowParticleNumber/step.h
        LowParticleNumber/tools.h
        )

find_package(GSL REQUIRED)
target_link_libraries(efficient_MDGFRD GSL::gsl GSL::gslcblas)