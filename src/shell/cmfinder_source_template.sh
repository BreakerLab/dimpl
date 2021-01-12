# Pull in variables and set up necessary executables
source scripts/cluster.conf

SEQDIR=$STEPNAME/$SEQNAME
if [ -f $SEQDIR/NO_HITS_FOUND ]; then echo "CMfinder did not need to run because prior CMsearch returned no hits."; exit; fi
BASECOMMAND="cmfinder04.pl -s1 7 -s2 7 -combine -fragmentary -skipClustalw -motifList $SEQDIR/motif.list --allCpus -commaSepEmFlags x--filter-non-frag,--max-degen-per-hit,2,--max-degen-flanking-nucs,7,--degen-keep,--amaa"
INITIALCOMMAND="$BASECOMMAND $SEQDIR/$SEQNAME.dedupe.fasta"
REALIGNCOMMAND="$($BASECOMMAND -justGetCmfinderCommand | cut -d':' -f2 | sed 's/\.$//' ) -o $SEQDIR/$SEQNAME.cmfinder.realign.sto -a $SEQDIR/$SEQNAME.cm.align.sto $SEQDIR/$SEQNAME.dedupe.fasta"
CMFINDERCOMMAND="$INITIALCOMMAND && $REALIGNCOMMAND"
