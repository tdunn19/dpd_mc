#!/bin/bash

set -o pipefail

dir='/Users/taylor/Documents/workspace/dpd_mc/data/v4.5/n'$1'd'$2'c'$3't'$4''


cp header.dat "$dir"/PQ
cd "$dir"/PQ

for u in {1..99}; do
    cat PQ"$u".dat >> header.dat
done

mv header.dat PQ.dat
mv PQ.dat ~/Documents/workspace/dpd_mc/analysis
cd ~/Documents/workspace/dpd_mc/analysis

set +e
./recon.exe < PQ.dat > F.dat

mv F.dat "$dir"
rm PQ.dat
echo 'test'
