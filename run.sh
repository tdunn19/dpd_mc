#!/bin/bash

nmon=10
density=1
job=1

dir=/home/tdunn

cd "$dir"/dpd/src/temp/'n'$nmon-'d'$density'-'$job''

make opt3
gcc -lm -o3 -o seed.exe seed.c 

echo -e 'n'$nmon'-d'$density'-'$job'' > output"$job".log

for u in {1..5}; do

	# Generate a random seed and add it to the input file
	sed -i -e "20d" dpd.inp
	./seed.exe >> dpd.inp

	echo 'Run '$u'.' >> output"$job".log
	./dpd.run >> output"$job".log

    mv energy.dat energy"$u".dat
    mv re.dat re"$u".dat
    mv rg.dat rg"$u".dat
    mv bond_length.dat bond_length"$u".dat
	
    mv energy"$u".dat "$dir"/scratch/dpd/n"$nmon"/d"$density"/energy
	mv re"$u".dat "$dir"/scratch/dpd/n"$nmon"/d"$density"/re
	mv rg"$u".dat "$dir"/scratch/dpd/n"$nmon"/d"$density"/rg
	mv bond_length"$u".dat "$dir"/scratch/dpd/n"$nmon"/d"$density"/bond_length

	echo 'Done '$u'.'
done

cp dpd.inp dpd"$job".inp 
mv dpd"$job".inp "$dir"/scratch/dpd/n"$nmon"/d"$density"/input

echo 'Done n'$nmon'-d'$density'-'$job'' >> output"$job".log
mv output"$job".log "$dir"/scratch/dpd/n"$nmon"/d"$density"/output
