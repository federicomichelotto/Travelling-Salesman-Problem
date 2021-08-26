#!/bin/bash

#SBATCH --job-name tsp_gen31
#SBATCH --output tsp_gen31_%j.out
#SBATCH --error tsp_gen31_%j.err

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

function run_instances() {

  for file in ../data/tests/meta/*; do
        echo "Submitted the run of GENETIC 3 for instance $file"
        srun ./tsp -f "$file" --meta -t 3600 --grasp 3 --popsize 1000 --offsize 100 --parsel 3 --sursel 1 -m 1 -v 2 -scores "gene_31"
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
