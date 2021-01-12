#!/bin/bash

thiscmd=`realpath $0`
parent_dir="$(dirname $thiscmd)"
batchfile=$parent_dir/scripts/blast_batchfile.sh
jobfile=$parent_dir/scripts/blast_jobfile.sh
mkdir -p $parent_dir/output

# Pull in variables and set up necessary executables
source $parent_dir/scripts/cluster.conf

echo "Generating batchfile at $batchfile"
dSQ.py --jobfile $jobfile --batch-file $batchfile --nice \
  --partition $PARTITION --mem 70G --chdir $parent_dir \
  --output output/blast_output%4a.out --job-name blast.$(basename $parent_dir)

sbatch $batchfile
