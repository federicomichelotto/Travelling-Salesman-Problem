#!/bin/bash

# Local variable
#SOLVER=("comp" "iter" "math" "heur" "tabu" "gene")
#GENETIC=("gen00" "gen01" "gen10" "gen11" "gen20" "gen21" "gen30" "gen31")
SOLVER=("comp")
path="/home/miliamikel/tsp/scripts"

for solver in "${SOLVER[@]}"; do

  if [ "$solver" == "comp" ]; then
    echo "Compact model submitted"
    for m in {1,2,5,8}; do
      sbatch "$path/$solver"32-tsp.slurm -m $m
    done

    for m in {3,4,6,7}; do
      sbatch "$path/$solver"64-tsp.slurm -m $m
    done

  elif [ "$solver" == "iter" ]; then
    echo "Iterative model submitted"

    for m in {9,10,11}; do
      sbatch "$path/$solver"-tsp.slurm -m $m
    done

  elif [ "$solver" == "math" ]; then
    echo "Matheuristic model submitted"

    for m in {0,1}; do
      sbatch "$path/$solver"-tsp.slurm -m $m
    done

  elif [ "$solver" == "heur" ]; then
    for m in {0,1,2,3}; do
      sbatch "$path/$solver"-tsp.slurm -m $m
    done

  elif [ "$solver" == "tabu" ]; then
    sbatch "$path/$solver"-tsp.slurm

  elif [ "$solver" == "gene" ]; then
    for genetic in "${GENETIC[@]}"; do
      sbatch "$path/$genetic"-tsp.slurm
    done

  else
    echo "Solver not found"

  fi
done
