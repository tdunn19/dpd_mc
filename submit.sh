#!/bin/bash

dir=/home/tdunn/scratch/dpd

# Remove parameter labels in input file
awk '{$1=""; print $0;}' run.inp > runtemp.inp

# Read input value and store parameters
density=$(sed -n '1p' runtemp.inp)
density=${density/ /} # Remove white space (for file names)
volume=$(sed -n '2p' runtemp.inp)
volume=${volume/ /}
nsteps=$(sed -n '3p' runtemp.inp)
nsteps=${nsteps/ /}
nfiles=$(sed -n '4p' runtemp.inp)
nfiles=${nfiles/ /}

rm runtemp.inp
cp run.sh "$dir"/src 
cp run.inp "$dir"/src 

cd "$dir"/src

# Change dpd.inp with parameters scanned from run.inp
sed -i 's/.*density.*/'$density'\t\t\tdensity/g' dpd.inp
sed -i 's/.*volume.*/'$volume'\t\t\tvolume/g' dpd.inp
sed -i 's/.*nsteps.*/'$nsteps'\t\t\tnsteps/g' dpd.inp


# Create directories if they don't already exist
mkdir "$dir"/results/
mkdir "$dir"/results/d"$density"
mkdir "$dir"/results/d"$density"/pressure
mkdir "$dir"/results/d"$density"/input
mkdir "$dir"/results/d"$density"/output

w=1
u=1

while [ $u -lt $nfiles ] 
do

	v=$((u+0)) # change

	loop='{'$u'..'$v'}; do'


	mkdir "$dir"/src/temp/'d'$density'-'$w''
	cd "$dir"/src
	cp * "$dir"/src/temp/'d'$density'-'$w''
	cp run.inp "$dir"/results/d"$density"/input
	
	cd "$dir"/src/temp/'d'$density'-'$w''
	sed -i "s/\(density *= *\).*/\1$density/" run.sh
	sed -i "s/\(part *= *\).*/\1$w/" run.sh
	sed -i "s/\(for u in *  *\).*/\1$loop/" run.sh

	qsub -l h_rt=1:00:00 run.sh # change

	echo 'd'$density'-'$w' from '$u' to '$v''
	w=`expr $w + 1`
	u=`expr $v + 1`
done

