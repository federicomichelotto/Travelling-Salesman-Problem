#!/bin/bash

cd cmake-build-debug || exit

for file in ../data/compact/*; do
  nodes=$(echo "$file" | grep -o -E '[0-9]+')
  if [ $nodes -le 14 ]; then
    for model in {0..5}; do
      if [ -f "$file" ]; then
        ./tsp -f "$file" -t 600 -m $model -v 2 -s $((1 + $RANDOM % 5000))
      fi
    done
  fi
done
