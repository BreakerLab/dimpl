#!/bin/bash

# Replace the line below with the path to the blast sequence database on your cluster 
export BLASTDB=/ysm-gpfs/datasets/db/blast
# If your cluster doesn't use `module load` replace with the necessary command to import Blast command line utils
module load BLAST+

#################################
# Do not modify below this line
#################################
BLASTCMD="blastx -evalue 0.01 -outfmt 5 -db nr -query"