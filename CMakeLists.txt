cmake_minimum_required(VERSION 3.9)
project(efficient_MDGFRD)

set(CMAKE_CXX_STANDARD 11)

add_executable(efficient_MDGFRD
        main.cpp
        LowParticleNumber/Fig4/Fig4.hpp
        LowParticleNumber/Fig9/Fig9.hpp
        LowParticleNumber/bruteForce.hpp
        LowParticleNumber/burst.hpp
        LowParticleNumber/checks.hpp
        LowParticleNumber/draw.hpp
        LowParticleNumber/greensFunct.hpp
        LowParticleNumber/init.hpp
        LowParticleNumber/parameters.hpp
        LowParticleNumber/print.hpp
        LowParticleNumber/run_aGF1.hpp
        LowParticleNumber/run_aGF2.hpp
        LowParticleNumber/run_BM.hpp
        LowParticleNumber/run_GF.hpp
        LowParticleNumber/run_hybrGF.hpp
        LowParticleNumber/shell.hpp
        LowParticleNumber/step.hpp
        LowParticleNumber/tools.hpp
        )

find_package(GSL REQUIRED)
target_link_libraries(efficient_MDGFRD GSL::gsl GSL::gslcblas)