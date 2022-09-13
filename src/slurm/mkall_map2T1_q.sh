#!/usr/bin/bash

declare -a lams=(0.1 1 3 10 30 100)
for lam in "${lams[@]}"; do
for tree in `seq 1 2`;do
./mkjob_map2T1.sh $lam $tree > "map_lam"$lam"tree"$tree".q"
sbatch "map_lam"$lam"tree"$tree".q"
done
done

