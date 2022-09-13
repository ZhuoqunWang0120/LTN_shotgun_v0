#!/usr/bin/bash

declare -a lams=(0.1 1 3 10 30 100)
for lam in "${lams[@]}"; do
./mkjob_tree_right.sh $lam > right$lam".q"
./mkjob_tree_left.sh $lam > left$lam".q"
./mkjob_tree_balanced.sh $lam > balanced$lam".q"
sbatch right$lam".q"
sbatch left$lam".q"
sbatch balanced$lam".q"
done

