#!/bin/bash
#  Glenn Gaffield - April 2020
#
#  Synchronize a sto file (by accession number) with a fasta file.
#    - Accession numbers present in the fasta file will be the only ones included in the sto.
#    - "#=GR" comment lines are removed, to clean up the output.
#    - If a fasta file is not provided, then only "#=GR" comment lines are removed.
#
#  Usage:  sync_and_clean_sto.sh {sto file} [fasta file] > {new sto file}

set -e

sto_file="$1"
fasta_file="$2"
if [ "$fasta_file" == "" ]; then
    fasta_file="/dev/null"
fi

( awk '/^>/ {print $1}' $fasta_file; cat $sto_file  ) |
    awk '/^>/ {sub("^>","",$1); accno[$1]=1; next} \
         /^#=GR/ {next};
         !/^[>#]/ && accno[$1] || "'$fasta_file'"=="/dev/null" {print; next} ;\
         /^[#/]/ || /^$/ {print}'
