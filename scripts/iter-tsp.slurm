#!/bin/bash

#SBATCH --job-name tsp_iter
#SBATCH --output tsp_iter_%A_%a.out
#SBATCH --error tsp_iter_%A_%a.err

#SBATCH --mail-type ALL
#SBATCH --mail-user mikele.milia@studenti.unipd.it

#SBATCH -n 1
#SBATCH -c 8
#SBATCH -p allgroups
#SBATCH -t 5-00:00:00
#SBATCH --mem 32G
#SBATCH --array=9,10,11

# Setup my personal project folder
DIRECTORY=/home/miliamikel/tsp
BUILD=build

function run_instances() {

  for file in ../data/tests/iter/*; do
    for i in {6..9}; do

        echo "Submitted the $i-th run of model ${SLURM_ARRAY_TASK_ID} for instance $file"
        srun ./tsp -f "$file" --ticks -t 1020000 -m "${SLURM_ARRAY_TASK_ID}" -v 4 -s $i -r $i

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
