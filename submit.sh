#!/bin/bash

dir=/home/tdunn/scratch/dpd

# Print out memory usage
du -h "$dir"/results
du -h "$dir"/src/temp
du -h ~/run.sh*
quota

echo "Abort submitting jobs? (y/n) "
read yn
if [ $yn == y ]; then
    echo 'Aborting..'
    exit
fi
echo 'Continuing..'

# Remove parameter labels in input file
awk '{$1=""; print $0;}' run.inp > runtemp.inp

# Read input value and store parameters
nmon=$(sed -n '1p' runtemp.inp)
nmon=${nmon/ /} # Remove white space (for file names)
density=$(sed -n '2p' runtemp.inp)
density=${density/ /}
volume=$(sed -n '3p' runtemp.inp)
volume=${volume/ /}
nsteps=$(sed -n '4p' runtemp.inp)
nsteps=${nsteps/ /}
nfiles=$(sed -n '5p' runtemp.inp)
nfiles=${nfiles/ /}

rm runtemp.inp
cp run.sh "$dir"/src 
cp run.inp "$dir"/src 

cd "$dir"/src

# Change dpd.inp and run.sh with parameters scanned from run.inp
sed -i "s/.*n_mon.*/$nmon\t\t\tn_mon/g" dpd.inp
sed -i "s/.*density.*/$density\t\t\tdensity/g" dpd.inp
sed -i "s/.*volume.*/$volume\t\t\tvolume/g" dpd.inp
sed -i "s/.*nsteps.*/$nsteps\t\t\tnsteps/g" dpd.inp
sed -i "s/\(nmon *= *\).*/\1$nmon/" run.sh
sed -i "s/\(density *= *\).*/\1$density/" run.sh

# Create directories if they don't already exist
mkdir "$dir"/results/
mkdir "$dir"/results/n"$nmon"
mkdir "$dir"/results/n"$nmon"/d"$density"
mkdir "$dir"/results/n"$nmon"/d"$density"/energy
mkdir "$dir"/results/n"$nmon"/d"$density"/re
mkdir "$dir"/results/n"$nmon"/d"$density"/rg
mkdir "$dir"/results/n"$nmon"/d"$density"/bond_length
mkdir "$dir"/results/n"$nmon"/d"$density"/input
mkdir "$dir"/results/n"$nmon"/d"$density"/output

w=1
u=1

while [ $u -lt $nfiles ] 
do

	v=$((u+1)) # change

	loop='{'$u'..'$v'}; do'

	mkdir "$dir"/src/temp/'n'$nmon'-d'$density'-'$w''
	cd "$dir"/src
	cp * "$dir"/src/temp/'n'$nmon'-d'$density'-'$w''
	
	cd "$dir"/src/temp/'n'$nmon'-d'$density'-'$w''
    sed -i "s/\(part *= *\).*/\1$w/" run.sh
	sed -i "s/\(for u in *  *\).*/\1$loop/" run.sh

	qsub -l h_rt=48:00:00 run.sh # change

	echo 'n'$nmon'-d'$density'-'$w' from '$u' to '$v''
	w=`expr $w + 1`
	u=`expr $v + 1`
done

cd "$dir"/src
cp run.inp "$dir"/results/n"$nmon"/d"$density"/input
rm run.sh run.inp
