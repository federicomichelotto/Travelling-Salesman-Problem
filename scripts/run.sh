#!/bin/bash

cd ../build || exit

for file in ../data/compact/*; do
  nodes=$(echo "$file" | grep -o -E '[0-9]+')
  if [ "$nodes" -le 110 ]; then # -eq , -le, -lt, -ge, -gt
    for model in {1..5}; do
      for i in {1..5}; do
        if [ -f "$file" ]; then
          ./tsp -f "$file" -t 600 -m $model -v 2 -s $i -r $i
        fi
      done
    done
  fi
done
