#!/bin/bash

# If your cluster doesn't use `module load`, add dSQ utils to path
module load dSQ
# Replace with the slurm partition your jobs will be run in
PARTITION=pi_breaker

#################################
# Do not modify below this line
#################################
thiscmd=`realpath $0`
parent_dir="$(dirname $thiscmd)"
mkdir -p $parent_dir/output
infernal_batchfile="$parent_dir/scripts/${STEPNAME}_infernal_batchfile.sh"
infernal_jobfile="$parent_dir/scripts/${STEPNAME}_infernal_jobfile.sh"

echo "Generating infernal batchfile at $infernal_batchfile"
dSQ.py --jobfile $infernal_jobfile --batch-file $infernal_batchfile -t 1-0\
  --partition $PARTITION --mem 8G --chdir $parent_dir --status-dir $parent_dir \
  --output output/infernal_output%4a.out --job-name infernal.$(basename $parent_dir)
infernal_RESPONSE=$(sbatch $infernal_batchfile)
infernal_JOBID=${infernal_RESPONSE##* }

# Fallback option to retry any jobs that fail due to error code 137 out of memory errors
infernal_bigmem_batchfile="$parent_dir/scripts/${STEPNAME}_infernal_bigmem_batchfile.sh"
infernal_bigmem_jobfile="$parent_dir/scripts/${STEPNAME}_infernal_bigmem_jobfile.sh"
infernal_bigmem_dsq="dSQ.py --jobfile $infernal_bigmem_jobfile --batch-file $infernal_bigmem_batchfile \
  --partition $PARTITION --mem 32G --chdir $parent_dir --status-dir $parent_dir -c 4 --time 2-0\
  --output output/infernal_bigmem_output%4a.out --job-name infernal_bigmem.$(basename $parent_dir) --submit"
bigmem_COMMAND=$'#!/bin/bash \n'"module load dSQ; cat $parent_dir/job_${infernal_JOBID}_status.tsv | cut -f 2,7 | grep ^137 | cut -f2 > $infernal_bigmem_jobfile; $infernal_bigmem_dsq"
bigmem_RESPONSE=$(echo "$bigmem_COMMAND" | sbatch -p $PARTITION --dependency afterany:${infernal_JOBID} -t 1-0 --job-name dsq_bigmem_infernal.$(basename $parent_dir))
bigmem_JOBID=${bigmem_RESPONSE##* }

# Fix coordinate discrepencies arising from searching a database of IGRs rather than genomes.
fix_COMMAND=$'#!/bin/bash\n find . -name \*.tblout | xargs -I {} sh -c "../scripts/fix_tblout_coords.sh {} > {}.tmp && mv {}.tmp {}"; find . -name \*.sto | xargs -I {} sh -c "../scripts/fix_sto_coords.sh {} > {}.tmp && mv {}.tmp {}"'
fix_RESPONSE=$(echo "$fix_COMMAND" | sbatch -p $PARTITION --dependency afterany:${bigmem_JOBID} -t 1-0 --job-name fixcoords.$(basename $parent_dir) --chdir $parent_dir/$STEPNAME --output ../output/fixcoords_${STEPNAME}.out)
fix_JOBID=${fix_RESPONSE##* }

# Convert all the results' .sto files into fasta files
fasta_COMMAND=$'#!/bin/bash\n module load HMMER; find . -type f -name \*.cm.align.sto -ls | awk \'($2>0){print $11}\' | xargs -I{} basename {} .cm.align.sto | xargs -I{} sh -c "esl-reformat fasta {}/{}.cm.align.sto > {}/{}.cmsearch.fasta"'
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
  --dependency afterany:${sample_JOBID} --cpus-per-task 2 --time 48:0:0

cmfinder_RESPONSE=$(sbatch $cmfinder_batchfile)
cmfinder_JOBID=${cmfinder_RESPONSE##* }
