#!/bin/bash
find . -type f -name \*.cmsearch.fasta -ls | \
xargs -I{} basename {} .cmsearch.fasta | \
xargs -I@@ bash -c "\
    awk 'BEGIN {RS=\">\";FS=\"\n\"} NR>1 {seq=\"\"; for (i=2;i<=NF;i++) seq=seq\$i; print \">\"\$1\"\t\"seq}' @@/@@.cmsearch.fasta |
    sort -u -k2,2 |
    tee >(shuf -n 1000 | awk 'BEGIN {FS=\"\t\"} {print \$1; print \$2}' > @@/@@.sample.fasta) |
    awk 'BEGIN {FS=\"\t\"} {print \$1; print \$2}' > @@/@@.dedupe.fasta"
