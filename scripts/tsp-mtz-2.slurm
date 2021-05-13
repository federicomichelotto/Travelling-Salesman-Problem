#!/bin/bash

#SBATCH --job-name tsp-mtz-2
#SBATCH --output output_%j.txt
#SBATCH --error errors_%j.txt

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

function run_instances() {
    echo "Submission : STARTING"

    for file in ../data/compact/*; do
      nodes=$(echo "$file" | grep -o -E '[0-9]+')
      if [[ $nodes -le 110 ]]; then # -eq , -le, -lt, -ge, -gt
        for model in {2..2}; do
          for i in {1..5}; do
            if [ -f "$file" ]; then
              echo "Submitted the $i run of model $model for instance $file"
              srun ./tsp -f "$file" --ticks -t 1020000 -m $model -v 4 -s $i -r $i
            fi
          done
        done
      fi
    done

    echo "Submission : DONE"
}

# Change the folder to my project directory
cd $DIRECTORY || exit

if [[ -d "$BUILD" ]]; then

  files=$(ls -A $BUILD)

  if [[ -z "$files" ]]; then

    echo "$BUILD directory is currently empty..."
    echo "Building executables : STARTING"

    cd $BUILD || exit

    cmake ..
    make
    echo "Building executables : DONE"

    run_instances

  else

    cd $BUILD || exit

    run_instances

  fi

else

  echo "Directory '$BUILD' does not exist, create it."

fi
