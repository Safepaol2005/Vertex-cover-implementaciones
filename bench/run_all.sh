#!/bin/bash
# Corre los 6 algoritmos (voraz, MIS R=1, MIS R=64, PD, DyV, HyVC) sobre los tres conjuntos, 2 procesos a la vez.
cd "$(dirname "$0")/.."
mk() { for f in $(ls $1/*.gr | sort); do
  echo "greedy $f 10"; echo "rmis $f 10 --runs 1 --seed 1"; echo "rmis $f 10 --runs 64 --seed 1"
  echo "dp $f 10"; echo "dc $f 10"; echo "hybrid $f 10"; done; }
mk data/pace2019/pace2019_track1_vc_exact_all > bench/jobs_pace2019.txt
mk data/synthetic > bench/jobs_synthetic.txt
mk data/pace2017 > bench/jobs_pace2017.txt
for ds in pace2019 synthetic pace2017; do
  xargs -P 2 -L 1 ./bench/run_one.sh < bench/jobs_$ds.txt > results/${ds}_raw.csv
done
