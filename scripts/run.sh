#!/bin/bash

# Local variable
#SOLVER=("comp" "iter" "math" "heur" "meta")
SOLVER=("comp")

# Move into build folder
for solver in "${SOLVER[@]}"; do

  if [ "$solver" == "comp" ]; then
    for m in {1..8}; do
      sbatch "$solver"-tsp.slurm -m $m
    done

  elif [ "$solver" == "iter" ]; then
    sbatch "$solver"-tsp.slurm

  elif [ "$solver" == "math" ]; then
    sbatch "$solver"-tsp.slurm

  elif [ "$solver" == "heur" ]; then
    sbatch "$solver"-tsp.slurm

  elif [ "$solver" == "meta" ]; then
    sbatch "$solver"-tsp.slurm

  else
    echo "Solver not found"

  fi
done
