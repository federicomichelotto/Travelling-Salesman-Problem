#!/bin/bash

# Local variable
SOLVER=("iter")
path="/home/miliamikel/tsp/scripts"

# Move into build folder
for solver in "${SOLVER[@]}"; do

  if [ "$solver" == "iter" ]; then
    sbatch "$path/$solver"-tsp-opt.slurm

  elif [ "$solver" == "heur" ]; then
    for m in {0,2,3}; do
      sbatch "$path/$solver"-tsp-opt.slurm -m "$m"
    done

  else
    echo "Solver not found"

  fi
done
