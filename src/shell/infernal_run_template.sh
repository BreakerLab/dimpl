#!/bin/bash


thiscmd=`realpath $0`
parent_dir="$(dirname $thiscmd)"
mkdir -p $parent_dir/output
infernal_batchfile="$parent_dir/scripts/${STEPNAME}_infernal_batchfile.sh"
infernal_jobfile="$parent_dir/scripts/${STEPNAME}_infernal_jobfile.sh"

# Pull in variables and set up necessary executables
source $parent_dir/scripts/cluster.conf

# Enable Email notification if EMAIL var defined in cluster.conf
if [ -n "$EMAIL" ]; then email_option="--mail-user $EMAIL --mail-type END"; fi

echo "Generating infernal batchfile at $infernal_batchfile"
dSQ.py --jobfile $infernal_jobfile --batch-file $infernal_batchfile \
  --partition $PARTITION --chdir $parent_dir --status-dir $parent_dir \
  --time 2-0 --cpus-per-task 7 --mem-per-cpu 8G \
  --output output/infernal_output%4a.out --job-name infernal.$(basename $parent_dir) > /dev/null
infernal_RESPONSE=$(sbatch $infernal_batchfile)
infernal_JOBID=${infernal_RESPONSE##* }

# Fix coordinate discrepencies arising from searching a database of IGRs rather than genomes.
fix_COMMAND=$'#!/bin/bash\n find . -name \*.tblout | xargs -I {} sh -c "../scripts/fix_tblout_coords.sh {} > {}.tmp && mv {}.tmp {}"; find . -name \*.sto | xargs -I {} sh -c "../scripts/fix_sto_coords.sh {} > {}.tmp && mv {}.tmp {}"; find . -name \*.tblout | while read tblout; do if [ `grep -cv "^#" $tblout` -eq 0 ]; then touch `dirname $tblout`/NO_HITS_FOUND; fi; done'
fix_RESPONSE=$(echo "$fix_COMMAND" | sbatch -p $PARTITION --dependency afterany:${infernal_JOBID} -t 1-0 --job-name fixcoords.$(basename $parent_dir) --chdir $parent_dir/$STEPNAME --output ../output/fixcoords_${STEPNAME}.out)
fix_JOBID=${fix_RESPONSE##* }

# Fallback option to retry any jobs that fail due to error code 137 out of memory errors
infernal_bigmem_batchfile="$parent_dir/scripts/${STEPNAME}_infernal_bigmem_batchfile.sh"
infernal_bigmem_jobfile="$parent_dir/scripts/${STEPNAME}_infernal_bigmem_jobfile.sh"
infernal_bigmem_dsq="set -x; bigmem_dsq_RESPONSE=\$(dSQ.py --jobfile $infernal_bigmem_jobfile --batch-file $infernal_bigmem_batchfile \
  --partition $PARTITION --mem 32G --chdir $parent_dir --status-dir $parent_dir -c 4 --time 2-0\
  --output output/infernal_bigmem_output%4a.out --job-name infernal_bigmem.$(basename $parent_dir) --submit | grep ^Submitted ); \
  bigmem_dsq_JOBID=\${bigmem_dsq_RESPONSE##* }; \
  scontrol update job ${fix_JOBID} Dependency=afterany:\${bigmem_dsq_JOBID}"
bigmem_COMMAND=$'#!/bin/bash \n'"module load dSQ; cat $parent_dir/job_${infernal_JOBID}_status.tsv | cut -f 2,7 | grep ^137 | cut -f2 > $infernal_bigmem_jobfile; if [ -s $infernal_bigmem_jobfile ]; then $infernal_bigmem_dsq; fi"
bigmem_RESPONSE=$(echo "$bigmem_COMMAND" | sbatch -p $PARTITION --dependency afterany:${infernal_JOBID} -t 1-0 --job-name dsq_bigmem_infernal.$(basename $parent_dir) --output output/dsq_bigmem_output.out)
bigmem_JOBID=${bigmem_RESPONSE##* }

scontrol update job ${fix_JOBID} Dependency=afterany:${bigmem_JOBID}

# Convert all the results' .sto files into fasta files
fasta_COMMAND=$'#!/bin/bash\n module load HMMER; sleep 20; find . -type f -name \*.cm.align.sto -ls | awk \'($2>0){print $11}\' | xargs -I{} basename {} .cm.align.sto | xargs -I{} sh -c "esl-reformat fasta {}/{}.cm.align.sto > {}/{}.cmsearch.fasta"'
fasta_RESPONSE=$(echo "$fasta_COMMAND" | sbatch -p $PARTITION --dependency afterany:${fix_JOBID} -t 1-0 --job-name sto2fasta.$(basename $parent_dir) --chdir $parent_dir/$STEPNAME --output ../output/sto2fasta_${STEPNAME}.out)
fasta_JOBID=${fasta_RESPONSE##* }

# For fasta files with more than 1000 hits, take a random sample of 1000.
sample_RESPONSE=$(sbatch -p $PARTITION --dependency afterany:${fasta_JOBID} -t 1-0 --job-name samplefasta.$(basename $parent_dir) --chdir $parent_dir/$STEPNAME --output ../output/samplefasta_${STEPNAME}.out $parent_dir/scripts/dedupe_sample.sh)
sample_JOBID=${sample_RESPONSE##* }


cmfinder_batchfile="$parent_dir/scripts/${STEPNAME}_cmfinder_batchfile.sh"
cmfinder_jobfile="$parent_dir/scripts/${STEPNAME}_cmfinder_jobfile.sh"

echo "Generating cmfinder batchfile at $cmfinder_batchfile"
dSQ.py --jobfile $cmfinder_jobfile --batch-file $cmfinder_batchfile \
  --partition $PARTITION --mem 16G --chdir $parent_dir --status-dir $parent_dir \
  --output output/cmfinder_output%4a.out --job-name cmfinder.$(basename $parent_dir) \
  --dependency afterany:${sample_JOBID} --cpus-per-task 2 --time 7-0 $email_option > /dev/null

cmfinder_RESPONSE=$(sbatch $cmfinder_batchfile)
cmfinder_JOBID=${cmfinder_RESPONSE##* }
