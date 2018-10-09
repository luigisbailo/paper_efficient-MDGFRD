for i in {0..10}
do
if [ -a $i.data ]
then
head -n +6 $i.data > temp.data
mv temp.data $i.data
fi
if [ -a $i.out ]
then
tail -n -6 $i.out > temp.out
mv temp.out $i.out
fi
done
