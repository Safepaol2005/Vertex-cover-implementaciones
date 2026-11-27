#!/bin/bash
# Genera las 35 instancias sintéticas del paper (semillas fijas).
cd "$(dirname "$0")/.." && mkdir -p data/synthetic && cd data/synthetic && G=../../bin/gen
for c in 1 2 3 4 6; do for s in 1 2 3; do $G er 1000 $c $s > er1000_c${c}_s$s.gr; done; done
for m in 2 3 5; do for s in 1 2; do $G ba 2000 $m $s > ba2000_m${m}_s$s.gr; done; done
for b in 8 12 16; do for s in 1 2; do $G band 1000 $b 0.4 $s > band1000_b${b}_s$s.gr; done; done
for k in 20 40 80 160; do $G trap $k > trap_k$k.gr; done
for c in 3 6; do for s in 1 2; do $G bip 500 500 $c $s > bip1000_c${c}_s$s.gr; done; done
