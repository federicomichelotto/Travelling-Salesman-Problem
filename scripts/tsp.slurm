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
SRC=src

# Change the folder to my project directory
cd $DIRECTORY || exit

# Change the folder to where the executable are located
if [[ -d "$SRC" ]]; then

  for file in $SRC/*; do
    echo "Checking $file ..."

    if [ "$(($(date +"%s") - $(stat -c "%Y" "$file")))" -gt "1800" ]; then
      echo "Updated recently file: $file"
      echo "Recompile needed"
      MAKE=1
      break
    fi
  done
fi

if [[ -d "$BUILD" ]]; then

  echo "Changing current directory to $BUILD ..."
  cd $BUILD || exit

  if [ $MAKE -eq 1 ]; then
    echo "Building the executables ..."
    cmake ..
    make
    echo "Building complete."
  fi

  echo "Starting the submission ..."

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
