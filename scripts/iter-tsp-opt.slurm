#!/bin/bash

#SBATCH --job-name tsp_iter_opt
#SBATCH --output tsp_iter_%j.out
#SBATCH --error tsp_iter_%j.err

#SBATCH --mail-type ALL
#SBATCH --mail-user mikele.milia@studenti.unipd.it

#SBATCH -n 1
#SBATCH -c 8
#SBATCH -p allgroups
#SBATCH -t 5-00:00:00
#SBATCH --mem 200G

# Setup my personal project folder
DIRECTORY=/home/miliamikel/tsp
BUILD=build

function run_instances() {

  for file in ../data/tests/iter/*; do
    for i in {6..9}; do

        echo "Submitted the $i-th run of model 11 for instance $file with 2-OPT"
        srun ./tsp -f "$file" --ticks --opt -t 1020000 -m 11 -v 4 -s "$(( i +10 ))" -r "$(( i +10 ))" -scores "iter_opt_$MODEL"

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
