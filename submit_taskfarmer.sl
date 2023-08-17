#!/bin/sh
#SBATCH -N 4 -c 128
#SBATCH -q regular
#SBATCH -t 05:00:00
#SBATCH -C cpu

export THREADS=128

runcommands.sh tasks.txt

