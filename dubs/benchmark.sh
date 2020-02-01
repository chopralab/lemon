#!/usr/bin/env bash
#PBS -d .
#PBS -l nodes=1:ppn=24
#PBS -l walltime=30:00

cp -a ../../full /dev/shm/

mkdir /tmp/dubs_out

module load gcc/6.3.0

for i in `seq 1 50`; do
    { time python3.6 dubs.py $test /dev/shm/full/ 1 /tmp/dubs_out/ > ${test}.out 2> ${test}.err ; } 2>> timings_2_$test
    rm /tmp/dubs_out/*
done

rm -fr /dev/shm/full
rm -fr /tmp/dubs_out

