#!/bin/bash

# Pull in variables and set up necessary executables
source scripts/cluster.conf

export BLASTDB
BLASTCMD="blastx -evalue 0.01 -outfmt 5 -db nr -query"
