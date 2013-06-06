#!/bin/bash

nmon=10
density=1
part=1

dir=/home/tdunn/scratch/dpd

cd "$dir"/src/temp/'n'$nmon-'d'$density'-'$part''

make opt3
gcc -lm -o3 -o seed.exe seed.c 

echo -e 'n'$nmon'-d'$density'-'$part'' > output"$part".log

for u in {1..100}; do

	# Generate a random seed and add it to the input file
	sed -i -e "20d" dpd.inp
	./seed.exe >> dpd.inp

	echo 'Run '$u'.' >> output"$part".log
	./dpd.run >> output"$part".log

    mv energy.dat energy"$u".dat
    mv re.dat re"$u".dat
    mv rg.dat rg"$u".dat
    mv bond_length.dat bond_length"$u".dat
	
	mv energy"$u".dat "$dir"/results/n"$nmon"/d"$density"/energy
	mv re"$u".dat "$dir"/results/n"$nmon"/d"$density"/re
	mv rg"$u".dat "$dir"/results/n"$nmon"/d"$density"/rg
	mv bond_length"$u".dat "$dir"/results/n"$nmon"/d"$density"/bond_length

	echo 'Done '$u'.'
done

cp dpd.inp dpd"$part".inp 
mv dpd"$part".inp "$dir"/results/n"$nmon"/d"$density"/input

echo 'Done n'$nmon'-d'$density'-'$part'' >> output"$part".log
mv output"$part".log "$dir"/results/n"$nmon"/d"$density"/output
