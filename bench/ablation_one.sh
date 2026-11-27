#!/bin/bash
# ablation_one.sh <flag> <archivo>  -> "flag,<línea de vcsolve>"
flag=$1; f=$2
echo "$flag,$("$(dirname "$0")/run_one.sh" hybrid "$f" 10 "$flag")"
