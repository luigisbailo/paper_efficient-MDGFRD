for i in {1..2}
do
head -n 5 perf$i.out > $i.data.txt
head -n 16 perf$i.out | tail -n 7 > $i.perf.txt
head -n 24 perf$i.out | tail -n 7 > $i.burst.txt
head -n 32 perf$i.out | tail -n 7 > $i.GF.txt
head -n 40 perf$i.out | tail -n 7 > $i.BM.txt
done

