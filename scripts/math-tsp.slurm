#!/bin/bash

#SBATCH --job-name tsp_math
#SBATCH --output tsp_math_%j.out
#SBATCH --error tsp_math_%j.err

#SBATCH --mail-type ALL
#SBATCH --mail-user mikele.milia@studenti.unipd.it

#SBATCH -n 1
#SBATCH -c 8
#SBATCH -p allgroups
#SBATCH -t 5-00:00:00
#SBATCH --mem 32G

# Setup my personal project folder
DIRECTORY=/home/miliamikel/tsp
BUILD=build

job_help() {
  # Display Help
  echo "Bash script that submit MAST simulation to SLURM scheduler"
  echo
  echo "Syntax: sbatch comp-tsp.slurm [-m|h]"
  echo "Options:"
  echo "-p     Set which model to use"
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

  for file in ../data/tests/math/*; do
    for i in {10..13}; do
        echo "Submitted the $i-th run of model $MODEL for instance $file"
        srun ./tsp -f "$file" --math --ticks -t 1020000 -m "$MODEL" -v 4 -s $i -r $i -scores "math_$MODEL"
    done
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
