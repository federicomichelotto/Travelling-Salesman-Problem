#!/bin/bash

# Local variable
SOLVER=("iter" "heur")

# Move into build folder
for solver in "${SOLVER[@]}"; do

  if [ "$solver" == "iter" ]; then
    sbatch "$solver"-tsp-opt.slurm

  elif [ "$solver" == "heur" ]; then
    for m in {0,1,2,3}; do
      sbatch "$solver"-tsp-opt.slurm -m $m
    done

  else
    echo "Solver not found"

  fi
done
