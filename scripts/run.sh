#!/bin/bash

# Local variable
#SOLVER=("comp" "iter" "math" "heur" "tabu" "gene")
GENETIC=("gen0" "gen1" "gen2" "gen3")
SOLVER=("comp")

# Move into build folder
for solver in "${SOLVER[@]}"; do

  if [ "$solver" == "comp" ]; then
    for m in {1,2,5,8}; do
      sbatch "$solver"-32-tsp.slurm -m $m
    done

    for m in {3,4,6,7}; do
      sbatch "$solver"-64-tsp.slurm -m $m
    done

  elif [ "$solver" == "iter" ]; then
    sbatch "$solver"-tsp.slurm

  elif [ "$solver" == "math" ]; then
    sbatch "$solver"-tsp.slurm

  elif [ "$solver" == "heur" ]; then
    sbatch "$solver"-tsp.slurm

  elif [ "$solver" == "tabu" ]; then
    sbatch "$solver"-tsp.slurm

  elif [ "$solver" == "gene" ]; then
    for genetic in "${GENETIC[@]}"; do
      sbatch "$genetic"-tsp.slurm
    done

  else
    echo "Solver not found"

  fi
done
