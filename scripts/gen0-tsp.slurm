#!/bin/bash

#SBATCH --job-name tsp_gen0
#SBATCH --output tsp_gen0_%j.out
#SBATCH --error tsp_gen0_%j.err

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
    for i in {101..102}; do
        echo "Submitted the $i-th run of GENETIC 0 for instance $file"
        srun ./tsp -f "$file" --meta -t 3600 --grasp 3 --popsize 1000 --offsize 100 --parsel 0 --sursel 0  -m 1 -v 4 -s $i -r $i
        srun ./tsp -f "$file" --meta -t 3600 --grasp 3 --popsize 1000 --offsize 100 --parsel 0 --sursel 1  -m 1 -v 4 -s "$(( i +10 ))" -r "$(( i +10 ))"
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
