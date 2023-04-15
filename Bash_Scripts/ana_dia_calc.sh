minmax ()
{

if [ -z "$1" ]
then
    echo "---parameter #1 is zero length---"
else
    echo "---parameter #1 is \"$1\"---"
fi

file=$1
fmin=(1 2 3);fmax=(1 2 3);fave=(1 2 3)
lo=(1 2 3);hi=(1 2 3);

ntot=$(sed -n '4,4p' $file | awk '{ print $1}')

for i in {6..8};
do
    n=$(($i - 6));
    lo[$n]=$(sed -n $i','$i'p' ${file} | awk '{ print $1}');
    hi[$n]=$(sed -n $i','$i'p' ${file} | awk '{ print $2}');
done

for i in {3..5}; 
do 

    n=$(($i - 3))   # array index start with 0

    fmin[$n]=$(sed '1,9d' ${file} | sort -gk${i} | head -n1 | awk -v c=$i '{ print $c}'); 
    fmax[$n]=$(sed '1,9d' ${file} | sort -gk${i} | tail -n1 | awk -v c=$i '{ print $c}'); 

    fave[$n]=$(sed '1,9d' ${file} | awk -v c=$i '{ sum += $c; n++} END {printf "%20.15e\n", sum/n}');

done



echo ${fmin[*]}
echo ${fmax[*]}

echo "output diameter"

for i in {0..2}
do
     echo ${fmin[$i]} ${fmax[$i]} | awk '{ print ($2-$1) }'

done
}

minmax $1
