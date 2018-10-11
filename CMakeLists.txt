cmake_minimum_required(VERSION 3.9)
project(Code)

set(CMAKE_CXX_STANDARD 11)

add_executable(Code
        2draft/alpha/fixAlpha/fixAlpha.hpp
        2draft/alpha/GF-BF/GF_BF.hpp
        2draft/diffusion/diffusion.hpp
        2draft/performance/performance.hpp
        2draft/bruteForce.hpp
        2draft/burst.hpp
        2draft/checks.hpp
        2draft/draw.hpp
        2draft/drawDef.hpp
        2draft/greensFunct.hpp
        2draft/greensFunctDef.hpp
        2draft/init.hpp
        2draft/main.cpp
        2draft/parameters.hpp
        2draft/print.hpp
        2draft/run_aGF1.hpp
        2draft/run_aGF2.hpp
        2draft/run_BM.hpp
        2draft/run_GF.hpp
        2draft/run_hybrGF.hpp
        2draft/shell.hpp
        2draft/step.hpp
        2draft/tools.hpp
#        2draftNeighb/performance/main.cpp
        2draftNeighb/bruteForce.hpp
        2draftNeighb/burst.hpp
        2draftNeighb/checks.hpp
        2draftNeighb/draw.hpp
        2draftNeighb/greensFunct.hpp
        2draftNeighb/init.hpp
#        2draftNeighb/main.cpp
        2draftNeighb/parameters.hpp
        2draftNeighb/print.hpp
        2draftNeighb/run_aGF.hpp
        2draftNeighb/run_BM.hpp
        2draftNeighb/run_GF.hpp
        2draftNeighb/shell.hpp
        2draftNeighb/step.hpp
        2draftNeighb/tools.hpp)

find_package(GSL REQUIRED)
target_link_libraries(Code GSL::gsl GSL::gslcblas)