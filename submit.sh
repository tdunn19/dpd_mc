#!/bin/bash

dir=/home/tdunn/dpd
scratch=/home/tdunn/scratch/dpd

# Print out memory usage
quota

echo "Abort submitting jobs? (y/n) "
read yn
if [ $yn == y ]; then
    echo 'Aborting..'
    exit
fi

du -h "$scratch"
echo "Delete scratch? (y/n) "
read yn
if [ $yn == y ]; then
    echo 'Deleting scratch..'
    rm -r "$scratch"
fi

du -h "$dir"/src/temp
echo "Delete temp? (y/n) "
read yn
if [ $yn == y ]; then
    echo 'Deleting temp folder..'
    rm -r "$dir"/src/temp
fi

mkdir "$dir"/src/temp
mkdir "$scratch"

# Remove parameter labels in input file
awk '{$1=""; print $0;}' run.inp > runtemp.inp

# Read input value and store parameters
nmon=$(sed -n '1p' runtemp.inp)
nmon=${nmon/ /} # Remove white space (for file names)
density_s=$(sed -n '2p' runtemp.inp)
density_s=${density_s/ /}
calc_list=$(sed -n '3p' runtemp.inp)
calc_list=${calc_list/ /}
length_x=$(sed -n '4p' runtemp.inp)
length_x=${length_x/ /}
length_y=$(sed -n '5p' runtemp.inp)
length_y=${length_y/ /}
length_z=$(sed -n '6p' runtemp.inp)
length_z=${length_z/ /}
n_steps=$(sed -n '7p' runtemp.inp)
n_steps=${n_steps/ /}
n_files=$(sed -n '8p' runtemp.inp)
n_files=${n_files/ /}
runs_per_job=$(sed -n '9p' runtemp.inp)
runs_per_job=${runs_per_job/ /}
time=$(sed -n '10p' runtemp.inp)
time=${time/ /}

if [ "$calc_list" = 1 ]; then
    calc="cell"
else
    calc="brute"
fi

rm runtemp.inp
cp run.sh run.inp "$dir"/src 

cd "$dir"/src
make opt3
gcc -lm -o3 -o seed.exe seed.c

# Change dpd.inp and run.sh with parameters scanned from run.inp
sed -i "s/.*n_mon.*/$nmon\t\t\tn_mon/g" dpd.inp
sed -i "s/.*density_s.*/$density_s\t\t\tdensity_s/g" dpd.inp
sed -i "s/.*calc_list.*/$calc_list\t\t\tcalc_list/g" dpd.inp
sed -i "s/.*length_x.*/$length_x\t\t\tlength_x/g" dpd.inp
sed -i "s/.*length_y.*/$length_y\t\t\tlength_y/g" dpd.inp
sed -i "s/.*length_z.*/$length_z\t\t\tlength_z/g" dpd.inp
sed -i "s/.*n_steps.*/$n_steps\t\t\tn_steps/g" dpd.inp

sed -i "s/\(nmon *= *\).*/\1$nmon/" run.sh
sed -i "s/\(density_s *= *\).*/\1$density_s/" run.sh
sed -i "s/\(calc *= *\).*/\1$calc/" run.sh

# Create directories to store data
mkdir -p "$scratch"/"$calc"/n"$nmon"/d"$density_s"
mkdir "$scratch"/"$calc"/n"$nmon"/d"$density_s"/energy
mkdir "$scratch"/"$calc"/n"$nmon"/d"$density_s"/re
mkdir "$scratch"/"$calc"/n"$nmon"/d"$density_s"/rg
mkdir "$scratch"/"$calc"/n"$nmon"/d"$density_s"/bond_length
mkdir "$scratch"/"$calc"/n"$nmon"/d"$density_s"/input
mkdir "$scratch"/"$calc"/n"$nmon"/d"$density_s"/output

w=1
u=1

while [ $u -lt $n_files ] 
do

	v=$((u+runs_per_job-1))

	loop='{'$u'..'$v'}; do'

	mkdir "$dir"/src/temp/'n'$nmon'-d'$density_s'-'$calc'-'$w''
	cd "$dir"/src
	cp dpd.run dpd.inp seed.exe run.sh "$dir"/src/temp/'n'$nmon'-d'$density_s'-'$calc'-'$w''
	
	cd "$dir"/src/temp/'n'$nmon'-d'$density_s'-'$calc'-'$w''
    sed -i "s/\(job *= *\).*/\1$w/" run.sh
	sed -i "s/\(for u in *  *\).*/\1$loop/" run.sh

	qsub -cwd -l h_rt="$time" run.sh

	echo 'n'$nmon'-d'$density_s'-'$calc'-'$w' from '$u' to '$v''
	w=`expr $w + 1`
	u=`expr $v + 1`
done

cd "$dir"/src
cp run.inp "$scratch"/"$calc"/n"$nmon"/d"$density_s"/input
rm run.sh run.inp dpd.run seed.exe
