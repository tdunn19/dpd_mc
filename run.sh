#!/bin/bash

nmon=10
density_s=1
calc="cell"
job=1

dir=/home/tdunn/dpd
scratch=/home/tdunn/scratch/dpd

cd "$dir"/src/temp/'n'$nmon-'d'$density_s'-'$calc'-'$job''

echo -e 'n'$nmon'-d'$density_s'-'$calc'-'$job'' > output"$job".log

for u in {1..5}; do

	# Generate a random seed and add it to the input file
	sed -i -e "30d" dpd.inp
	./seed.exe >> dpd.inp

	echo 'Run '$u'.' >> output"$job".log
	./dpd.run >> output"$job".log

    mv energy.dat energy"$u".dat
    mv re.dat re"$u".dat
    mv rg.dat rg"$u".dat
    mv bond_length.dat bond_length"$u".dat
	
    mv energy"$u".dat "$scratch"/"$calc"/n"$nmon"/d"$density_s"/energy
	mv re"$u".dat "$scratch"/"$calc"/n"$nmon"/d"$density_s"/re
	mv rg"$u".dat "$scratch"/"$calc"/n"$nmon"/d"$density_s"/rg
	mv bond_length"$u".dat "$scratch"/"$calc"/n"$nmon"/d"$density_s"/bond_length

	echo 'Done '$u'.'
done

cp dpd.inp dpd"$job".inp 
mv dpd"$job".inp "$scratch"/"$calc"/n"$nmon"/d"$density_s"/input

echo 'Done n'$nmon'-d'$density_s'-'$calc'-'$job'' >> output"$job".log
mv output"$job".log "$scratch"/"$calc"/n"$nmon"/d"$density_s"/output
