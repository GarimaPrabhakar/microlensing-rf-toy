#!/usr/bin/env bash

module load python
module load conda

cd "/pscratch/sd/g/garimap/microlensingToy/"

python generate_training_set.py $1 $2 $3
