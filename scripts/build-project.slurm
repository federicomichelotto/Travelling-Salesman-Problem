#!/bin/bash

#SBATCH --job-name cmake
#SBATCH --output tsp_cmake_%j.out
#SBATCH --error tsp_cmake_%j.err

#SBATCH --mail-type ALL
#SBATCH --mail-user mikele.milia@studenti.unipd.it

#SBATCH -n 1
#SBATCH -c 8
#SBATCH -p allgroups
#SBATCH -t 00:05:00
#SBATCH --mem 1G

# Setup my personal project folder
DIRECTORY=/home/miliamikel/tsp
BUILD=build

# Change the folder to my project directory
cd $DIRECTORY || exit

if [[ -d "$BUILD" ]]; then

  cd $BUILD || exit

  echo "$BUILD directory is currently empty..."
  echo "Building executables..."

  cmake ..
  make

else

  echo "Directory '$BUILD' does not exist, create it."

fi
