#!/bin/bash

dir=/home/tdunn

# Print out memory usage
quota

echo "Abort submitting jobs? (y/n) "
read yn
if [ $yn == y ]; then
    echo 'Aborting..'
    exit
fi

du -h "$dir"/scratch
echo "Delete scratch? (y/n) "
read yn
if [ $yn == y ]; then
    echo 'Deleting scratch..'
    rm -r "$dir"/scratch
fi

du -h "$dir"/dpd/src/temp
echo "Delete temp? (y/n) "
read yn
if [ $yn == y ]; then
    echo 'Deleting temp folder..'
    rm -r "$dir"/dpd/src/temp
fi

mkdir "$dir"/dpd/src/temp
mkdir "$dir"/scratch/dpd

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
runs_per_job=$(sed -n '6p' runtemp.inp)
runs_per_job=${runs_per_job/ /}
time=$(sed -n '7p' runtemp.inp)
time=${time/ /}

rm runtemp.inp
cp run.sh "$dir"/dpd/src 
cp run.inp "$dir"/dpd/src 

cd "$dir"/dpd/src

# Change dpd.inp and run.sh with parameters scanned from run.inp
sed -i "s/.*n_mon.*/$nmon\t\t\tn_mon/g" dpd.inp
sed -i "s/.*density.*/$density\t\t\tdensity/g" dpd.inp
sed -i "s/.*volume.*/$volume\t\t\tvolume/g" dpd.inp
sed -i "s/.*nsteps.*/$nsteps\t\t\tnsteps/g" dpd.inp
sed -i "s/\(nmon *= *\).*/\1$nmon/" run.sh
sed -i "s/\(density *= *\).*/\1$density/" run.sh
sed -i "s/\(runs_per_job *= *\).*/\1$runs_per_job/" run.sh

# Create directories if they don't already exist
mkdir "$dir"/scratch/dpd/n"$nmon"
mkdir "$dir"/scratch/dpd/n"$nmon"/d"$density"
mkdir "$dir"/scratch/dpd/n"$nmon"/d"$density"/energy
mkdir "$dir"/scratch/dpd/n"$nmon"/d"$density"/re
mkdir "$dir"/scratch/dpd/n"$nmon"/d"$density"/rg
mkdir "$dir"/scratch/dpd/n"$nmon"/d"$density"/bond_length
mkdir "$dir"/scratch/dpd/n"$nmon"/d"$density"/input
mkdir "$dir"/scratch/dpd/n"$nmon"/d"$density"/output

w=1
u=1

while [ $u -lt $nfiles ] 
do

	v=$((u+runs_per_job-1))

	loop='{'$u'..'$v'}; do'

	mkdir "$dir"/dpd/src/temp/'n'$nmon'-d'$density'-'$w''
	cd "$dir"/dpd/src
	cp * "$dir"/dpd/src/temp/'n'$nmon'-d'$density'-'$w''
	
	cd "$dir"/dpd/src/temp/'n'$nmon'-d'$density'-'$w''
    sed -i "s/\(job *= *\).*/\1$w/" run.sh
	sed -i "s/\(for u in *  *\).*/\1$loop/" run.sh

	qsub -cwd h_rt="$time" run.sh

	echo 'n'$nmon'-d'$density'-'$w' from '$u' to '$v''
	w=`expr $w + 1`
	u=`expr $v + 1`
done

cd "$dir"/dpd/src
cp run.inp "$dir"/results/n"$nmon"/d"$density"/input
rm run.sh run.inp
