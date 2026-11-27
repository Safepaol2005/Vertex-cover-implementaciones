#!/bin/bash
# run_one.sh <alg> <archivo> <limite_s> [args extra]  -> una línea CSV (o una línea de error)
alg=$1; f=$2; T=$3; shift 3
out=$(timeout $((${T%.*} + 15)) "$(dirname "$0")/../bin/vcsolve" "$alg" "$f" --time "$T" "$@" 2>/dev/null)
rc=$?
if [ -z "$out" ]; then echo "$alg,$f,,,,,killed_rc$rc,"; else echo "$out"; fi
