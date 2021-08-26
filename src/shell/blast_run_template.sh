#!/bin/bash

thiscmd=`realpath $0`
parent_dir="$(dirname $thiscmd)"
batchfile=$parent_dir/scripts/blast_batchfile.sh
jobfile=$parent_dir/scripts/blast_jobfile.sh
mkdir -p $parent_dir/output

# Pull in variables and set up necessary executables
source $parent_dir/scripts/cluster.conf

# Enable Email notification if EMAIL var defined in cluster.conf
if [ -n "$EMAIL" ]; then email_option="--mail-user $EMAIL --mail-type END"; fi

echo "Generating batchfile at $batchfile"
dSQ.py --jobfile $jobfile --batch-file $batchfile --nice \
  --partition $PARTITION --mem 70G --chdir $parent_dir \
  --output output/blast_output%4a.out --job-name blast.$(basename $parent_dir) $email_option

sbatch $batchfile
