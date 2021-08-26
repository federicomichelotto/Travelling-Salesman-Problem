#!/bin/bash

#SBATCH --job-name tsp_heur
#SBATCH --output tsp_heur_%j.out
#SBATCH --error tsp_heur_%j.err

#SBATCH --mail-type ALL
#SBATCH --mail-user federico.michelotto@studenti.unipd.it

#SBATCH -n 1
#SBATCH -c 8
#SBATCH -p allgroups
#SBATCH -t 5-00:00:00
#SBATCH --mem 16GB

# Setup my personal project folder
DIRECTORY=/home/michelfede/OR2/travelling-salesman-problem
BUILD=build

job_help() {
  # Display Help
  echo "Bash script that submit MAST simulation to SLURM scheduler"
  echo
  echo "Syntax: sbatch comp-tsp.slurm [-m|h]"
  echo "Options:"
  echo "-m     Set which model to use"
  echo "-h     Print the help."
}

# Retrieve custom project path
while getopts m:h option; do
  case "${option}" in
  m) MODEL=${OPTARG} ;;
  h) job_help
     exit ;;
  *) echo "Command not supported"
     exit ;;
  esac
done

function run_instances() {

  for file in ../data/tests/heur/*; do
        echo "Submitted the run of model $MODEL for instance $file"
        if [ "$MODEL" == 0 ]; then # NN
          srun ./tsp -f "$file" --heur -t 3600 -m "$MODEL" -scores "heur_0"
          srun ./tsp -f "$file" --heur --grasp 3 -t 3600 -m "$MODEL" -scores "heur_0_grasp"
        fi

        if [ "$MODEL" == 1 ]; then #NN + GRASP + 2-OPT
          srun ./tsp -f "$file" --heur --grasp 3 -t 3600 -m "$MODEL" -scores "heur_1"
        fi

	if [ "$MODEL" == 3 ]; then # EXTRA MILEAGE with furthest starting nodes
          srun ./tsp -f "$file" --heur -t 3600 -m "$MODEL" -scores "heur_3"
        fi
  done

}

# Change the folder to my project directory
cd $DIRECTORY || exit

if [[ -d "$BUILD" ]]; then

  files=$(ls -A $BUILD)

  cd $BUILD || exit

  if [[ -z "$files" ]]; then

    echo "$BUILD directory is currently empty..."
    echo "Building executables..."

    cmake ..
    make

  fi

  # Running instances
  run_instances

else

  echo "Directory '$BUILD' does not exist, create it."

fi
