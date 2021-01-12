# Pull in variables and set up necessary executables
source scripts/cluster.conf

SEQDIR=$STEPNAME/$SEQNAME
CMBUILDCOMMAND="cmbuild -F $NOSS -o $SEQDIR/$SEQNAME.out  $SEQDIR/$SEQNAME.cm  $SEQDIR/$SEQNAME.sto"
CMCALIBRATECOMMAND="cmcalibrate  $SEQDIR/$SEQNAME.cm"
CMSEARCHCOMMAND="cmsearch -E 0.01 -o  $SEQDIR/$SEQNAME.cm.out --tblout  $SEQDIR/$SEQNAME.cm.tblout -A  $SEQDIR/$SEQNAME.cm.align.sto  $SEQDIR/$SEQNAME.cm $DATABASE"
