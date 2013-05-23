#!/bin/bash

density=1
part=1

dir=/home/tdunn/scratch/dpd

cd "$dir"/src/temp/'d'$density'-'$part''

# Remove parameter labels in input file
awk '{$1=""; print $0;}' run.inp > runtemp.inp

# Read input value and store parameters
density=$(sed -n '1p' runtemp.inp)
density=${density/ /}
volume=$(sed -n '2p' runtemp.inp)
volume=${volume/ /}
nsteps=$(sed -n '3p' runtemp.inp)
nsteps=${nsteps/ /}
nfiles=$(sed -n '4p' runtemp.inp)
nfiles=${nfiles/ /}

rm runtemp.inp

make opt3
gcc -lm -o3 -o seed.exe seed.c 

echo -e 'd'$density'-'$part'' > output"$part".log

for u in {1..100}; do

	# Generate a random seed and add it to the input file
	sed -i -e "9d" dpd.inp
	./seed.exe >> dpd.inp

	echo 'Run '$u'.' >> output"$part".log
	./dpd.run >> output"$part".log
	
	mv pressure.dat pressure"$u".dat
	
	mv pressure"$u".dat "$dir"/results/d"$density"/pressure

	echo 'Done '$u'.'
done

cp dpd.inp dpd"$part".inp 
mv dpd"$part".inp "$dir"/results/d"$density"/input

echo 'Done d'$density'-'$part'' >> output"$part".log
mv output"$part".log "$dir"/results/d"$density"/output
