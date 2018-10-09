./Perf/main 0.01 2.5 0.01 2.5 0.1 5 > Perf/Perf1.out
echo 1Done
./Perf/main 0.01 1.5 0.001 3.5 0.1 5 > Perf/Perf2.out 
echo 2Done
./Comp/main 0.01 2.5 0.01 2.5 0.1 5 > Comp/Perf.out
echo 3Done
