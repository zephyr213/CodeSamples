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

}

declare -a f1min f1max f1ave f2min f2max f2ave lo1 hi1 lo2 hi2
# processing str1

## input files ##########

infile1=~/data/npgen3/npeq_3511p3/traj/current.atom.210000 
infile2=~/data/npgen3/npeq_3511p3/traj/current.atom.210000

infile1=./i8391v8391p3
infile2=./i8391v8391p3

# -- input paramters
d0=5    # initialdistance between two nps
dim=1   # join dimension, (0 1 2) corresponds to (x y z)

########################

minmax $infile1
for i in {0..2}
do
    f1min[$i]=${fmin[$i]}
    f1max[$i]=${fmax[$i]}
    f1ave[$i]=${fave[$i]}
    lo1[$i]=${lo[$i]}
    hi1[$i]=${hi[$i]}
    nt1=$ntot
done
# processing str2

minmax $infile2
for i in {0..2}
do
    f2min[$i]=${fmin[$i]}
    f2max[$i]=${fmax[$i]}
    f2ave[$i]=${fave[$i]}
    lo2[$i]=${lo[$i]}
    hi2[$i]=${hi[$i]}
    nt2=$ntot
done

echo "np1 centeroid:"
echo ${f1ave[*]}

echo "np2 centeroid:"
echo ${f2ave[*]}

# define a few variables

dimindex=(x y z)

echo "---Begin join---"
d1=$(echo ${f2max[$dim]} ${f1min[$dim]} ${d0} | awk '{ printf "%20.15e", ($1-$2+$3) }')
echo "d1:"$d1

displ=(0 0 0)
displ[$dim]=$d1
dx=${displ[0]};dy=${displ[1]};dz=${displ[2]}


ofile=str.out
str1=$(echo "$infile1" | grep -o -P '.{0,6}p3')
str2=$(echo "$infile2" | grep -o -P '.{0,6}p3')

nt=$(echo $nt1 $nt2 | awk '{ print ($1+$2) }')

echo $nt
echo "---Output now---"

echo $str1 $str2 > $ofile
echo >> $ofile
echo $nt " atoms" >> $ofile
echo "1 atom types" >> $ofile

for i in {0..2}
do
    dlo1=${lo1[$i]};dlo2=${lo2[$i]};
    dhi1=${hi1[$i]};dhi2=${hi2[$i]};
    if [ "$i" -ne "$dim" ]
    then
	lo0=$(awk -v a="$dlo1" -v b="$dlo2" 'BEGIN{print (a>b)?b:a}')
	hi0=$(awk -v a="$dhi1" -v b="$dhi2" 'BEGIN{print (a<b)?b:a}')
    else
	lo0=$(echo $dlo2 $d1 | awk '{ printf "%20.15e", ($1-$2) }')
	hi0=$dhi1
    fi
    echo $lo0 $hi0 ${dimindex[$i]}lo ${dimindex[$i]}hi >> $ofile
done
echo >> $ofile
echo "Atoms" >> $ofile
echo >> $ofile

sed '1,9d' $infile1 >> $ofile
sed '1,9d' $infile2 | awk -v d1=$dx -v d2=$dy -v d3=$dz -v nl=$nt1 '{ printf "%10d  %2d  %-1.16e %-1.16e %-1.16e\n", ($1+nl),($2),($3-d1),($4-d2),($5-d3) }' >> $ofile
