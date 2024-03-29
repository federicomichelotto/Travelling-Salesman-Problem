#!/bin/bash

#SBATCH --job-name tsp_tabu
#SBATCH --output tsp_tabu_%j.out
#SBATCH --error tsp_tabu_%j.err

#SBATCH --mail-type ALL
#SBATCH --mail-user mikele.milia@studenti.unipd.it

#SBATCH -n 1
#SBATCH -c 8
#SBATCH -p allgroups
#SBATCH -t 5-00:00:00
#SBATCH --mem 16GB

# Setup my personal project folder
DIRECTORY=/home/miliamikel/tsp
BUILD=build

function run_instances() {

  for file in ../data/tests/meta/*; do
        echo "Submitted the run of TABU SEARCH for instance $file"
	srun ./tsp -f "$file" --meta -t 3600 -m 0 -s "$i" -scores "tabu"
	srun ./tsp -f "$file" --meta -t 3600 --grasp 3 -m 0 -s "$i" -scores "tabu_grasp"
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
