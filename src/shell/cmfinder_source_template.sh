# Replace path below with the location of the cmfinder executables on your compute cluster
export PATH=$PATH:/ysm-gpfs/pi/breaker/software/packages/cmfinder-0.4.1.18/bin

#################################
# Do not modify below this line
#################################
SEQDIR=$STEPNAME/$SEQNAME
BASECOMMAND="cmfinder04.pl -s1 7 -s2 7 -combine -fragmentary -skipClustalw -motifList $SEQDIR/motif.list --allCpus -commaSepEmFlags x--filter-non-frag,--max-degen-per-hit,2,--max-degen-flanking-nucs,7,--degen-keep,--amaa -commaSepSummarizeFlags x--fragmentary,--degen-keep"
INITIALCOMMAND="$BASECOMMAND $SEQDIR/$SEQNAME.dedupe.fasta"
REALIGNCOMMAND="$($BASECOMMAND -justGetCmfinderCommand | cut -d':' -f2 | sed 's/\.$//' ) -o $SEQDIR/$SEQNAME.cmfinder.realign.sto -a $SEQDIR/$SEQNAME.cm.align.sto $SEQDIR/$SEQNAME.dedupe.fasta"
CMFINDERCOMMAND="$INITIALCOMMAND; $REALIGNCOMMAND"