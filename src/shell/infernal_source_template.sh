# If your cluster doesn't use `module load` replace with the necessary command to import Infernal command line utils
module load Infernal

# Replace with the path to the downloaded microbial igr database
DATABASE=/gpfs/ysm/pi/breaker/data/refseq98/s50.igr.fasta

#################################
# Do not modify below this line
#################################
SEQDIR=$STEPNAME/$SEQNAME
CMBUILDCOMMAND="cmbuild -F $NOSS -o $SEQDIR/$SEQNAME.out  $SEQDIR/$SEQNAME.cm  $SEQDIR/$SEQNAME.sto"
CMCALIBRATECOMMAND="cmcalibrate  $SEQDIR/$SEQNAME.cm"
CMSEARCHCOMMAND="cmsearch -E 0.01 -o  $SEQDIR/$SEQNAME.cm.out --tblout  $SEQDIR/$SEQNAME.cm.tblout -A  $SEQDIR/$SEQNAME.cm.align.sto  $SEQDIR/$SEQNAME.cm $DATABASE"