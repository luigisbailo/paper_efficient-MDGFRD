for i in {1..2}
do
head -n 6 perf$i.out > $i.data.txt
head -n 14 perf$i.out | tail -n 4 > $i.perf.txt
head -n 19 perf$i.out | tail -n 4 > $i.burst.txt
head -n 24 perf$i.out | tail -n 4 > $i.GF.txt
head -n 29 perf$i.out | tail -n 4 > $i.BM.txt
done

