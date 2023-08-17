#!/bin/sh
#SBATCH -N 4 -c 128
#SBATCH -q debug
#SBATCH -t 00:30:00
#SBATCH -C cpu

export THREADS=128

runcommands.sh tasks.txt

