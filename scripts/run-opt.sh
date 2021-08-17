#!/bin/bash

# Local variable
SOLVER=("iter" "heur")

# Move into build folder
for solver in "${SOLVER[@]}"; do

  if [ "$solver" == "iter" ]; then
    sbatch "$solver"-tsp-opt.slurm

  elif [ "$solver" == "heur" ]; then
    sbatch "$solver"-tsp-opt.slurm

  else
    echo "Solver not found"

  fi
done
