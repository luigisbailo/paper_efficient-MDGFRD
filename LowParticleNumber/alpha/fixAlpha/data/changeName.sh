for i in {0..9}
do
if [ -a $i.out ]
then
mv $i.out $i.out
cp $i.out $i.data
fi
done

