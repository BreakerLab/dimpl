#!/bin/bash

# If your cluster doesn't use `module load`, add dSQ utils to path
module load dSQ
# Replace with the slurm partition your jobs will be run in
PARTITION=pi_breaker

#################################
# Do not modify below this line
#################################
thiscmd=`realpath $0`
parent_dir="$(dirname $thiscmd)"
batchfile=$parent_dir/scripts/blast_batchfile.sh
jobfile=$parent_dir/scripts/blast_jobfile.sh
mkdir -p $parent_dir/output

echo "Generating batchfile at $batchfile"
dSQ.py --jobfile $jobfile --batch-file $batchfile --nice \
  --partition $PARTITION --mem 70G --chdir $parent_dir \
  --output output/blast_output%4a.out --job-name blast.$(basename $parent_dir)

sbatch $batchfile
