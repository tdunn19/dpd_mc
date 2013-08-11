#!/bin/bash


# Directories
dir=/home/tdunn/dpd_mc
scratch=/home/tdunn/scratch/dpd_mc

# Print out memory usage
quota
echo "Abort submitting jobs? (y/n) "
read yn
if [ $yn == y ]; then
    echo 'Aborting..'
    exit
fi

# Cleanup scratch folder
du -h "$scratch"
echo "Delete scratch? (y/n) "
read yn
if [ $yn == y ]; then
    echo 'Deleting scratch..'
    rm -r "$scratch"
fi

# Cleanup temp folder
du -h "$dir"/src/temp
echo "Delete temp? (y/n) "
read yn
if [ $yn == y ]; then
    echo 'Deleting temp folder..'
    rm -r "$dir"/src/temp
fi

# Exit early if an error code is returned
set -e
set -o pipefail

# Read run.inp
n_mon=$(awk '/n_mon/ {print $2}' run.inp)
density_s=$(awk '/density_s/ {print $2}' run.inp)
n_wins=$(awk '/n_wins/ {print $2}' run.inp)
a_ms_cis=$(awk '/a_ms_cis/ {print $2}' run.inp)
a_ms_trans=$(awk '/a_ms_trans/ {print $2}' run.inp)
n_steps=$(awk '/n_steps/ {print $2}' run.inp)
calc_list=$(awk '/calc_list/ {print $2}' run.inp)
runs_per_job=$(awk '/runs_per_job/ {print $2}' run.inp)
run_time=$(awk '/run_time/ {print $2}' run.inp)

# Modify dpd.inp
sed -i "s/.*n_mon.*/$n_mon\t\t\tn_mon/g" src/dpd.inp
sed -i "s/.*density_s.*/$density_s\t\t\tdensity_s/g" src/dpd.inp
sed -i "s/.*n_wins.*/$n_wins\t\t\tn_wins/g" src/dpd.inp
sed -i "s/.*a_ms_cis.*/$a_ms_cis\t\t\ta_ms_cis/g" src/dpd.inp
sed -i "s/.*a_ms_trans*/$a_ms_trans\t\t\ta_ms_trans/g" src/dpd.inp
sed -i "s/.*n_steps.*/$n_steps\t\t\tn_steps/g" src/dpd.inp
sed -i "s/.*calc_list.*/$calc_list\t\t\tcalc_list/g" src/dpd.inp

# Create directories
temp="$dir"/src/temp/n"$n_mon"d"$density_s"c"$a_ms_cis"t"$a_ms_trans"
data="$scratch"/n"$n_mon"d"$density_s"c"$a_ms_cis"t"$a_ms_trans"
mkdir -p "$temp"
mkdir -p "$data"
mkdir "$data"/PQ
mkdir "$data"/Q
mkdir "$data"/energy
mkdir "$data"/re
mkdir "$data"/re/cis
mkdir "$data"/re/trans
mkdir "$data"/rg
mkdir "$data"/rg/cis
mkdir "$data"/rg/trans
mkdir "$data"/input
mkdir "$data"/output

# Compile
make opt3 -C src
gcc -lm -o3 -o src/seed.exe src/seed.c

# Submit the jobs
for ((window=1,job=1;window<=$n_wins;job++)); do

    mkdir "$temp"/j"$job"
    cp src/dpd.run src/dpd.inp src/seed.exe run.sh "$temp"/j"$job"

    # Short.q <= 48:00:00
    # Medium.q <= 336:00:00
    # Long.q <= 720:00:00
    cd "$temp"/j"$job"
    qsub -cwd -l h_vmem=1G -l h_rt="$run_time" run.sh $job $window $((window+runs_per_job)) $data
    cd "$dir"

    echo 'Submitting n'$n_mon'd'$density_s'c'$a_ms_cis't'$a_ms_trans'j'$job' (w'$window' to w'$((window+runs_per_job-1))') for '$run_time''

    window=$((window+runs_per_job))
done

cp run.inp "$data"/input
rm src/dpd.run src/seed.exe
