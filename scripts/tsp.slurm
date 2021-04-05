#!/bin/bash

#SBATCH --job-name tsp
#SBATCH --output output.txt
#SBATCH --error errors.txt

#SBATCH --mail-type ALL
#SBATCH --mail-user mikele.milia@studenti.unipd.it

#SBATCH -n 1
#SBATCH -c 8
#SBATCH -p allgroups
#SBATCH -t 3-00:00:00
#SBATCH --mem 16G

# Setup my personal project folder
DIRECTORY=/home/miliamikel/tsp
BUILD=build

# Change the folder to my project directory
cd $DIRECTORY || exit

# Change the folder to where the executable are located

if [[ -d "$BUILD" ]]; then

  cd $BUILD || exit
  echo "Directory '$BUILD' exists"

  #  if [ "$(ls -A $BUILD)" ]; then
  echo "Building the executables ..."
  cmake ..
  make
  echo "Building complete."
  #  else
  echo "Starting the submission ..."
  #  fi

  for file in ../data/compact/*; do
    nodes=$(echo "$file" | grep -o -E '[0-9]+')
    if [[ $nodes -le 110 ]]; then # -eq , -le, -lt, -ge, -gt
      for model in {1..5}; do
        for i in {1..5}; do
          if [ -f "$file" ]; then
            echo "Submitted the $i run of model $model for instance $file"
            srun ./tsp -f "$file" -t 600 -m $model -v 2 -s $i -r $i
          fi
        done
      done
    fi
  done

  echo "Ending the submission ..."

else

  echo "Directory '$BUILD' does not exist, create it."

fi
