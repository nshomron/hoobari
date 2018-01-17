#!/bin/bash

## batch-queuing system arguments
#$ -N jobname
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -o .$JOB_NAME.$JOB_ID.log

# terminate on error
set -e

# region is given as an argument as described in freebayes' or hoobari's documentation
region=$1
out_dir=/out/directory
parents_vcf=parents.vcf.gz
cfdna_bam=cfdna.bam

output_vcf=${out_dir}/${region}.vcf
mkdir -p ${out_dir}
mkdir -p tmp_hb

# 1. Don't forget to load required modules for freebayes and python
# 2. Replace /path/to/ with the correct location, or delete it if the location is in your system's path.
# 3. 2>&1 redirects stderr to stdout, the the merged stream is split again: stdout is written to the
# output vcf and the remaining part (which is originally from stderr), is piped to hoobari's patch that
# creates the databases (and prints stderr for debugging).
# 4. To understand hoobari's patch parameters, run: python freebayes_hoobari_patch.py -h

/path/to/freebayes \
-d \
--region $region \
--fasta-reference ~/nshomron_tomr/tools/GATK/bundle/ucsc/ucsc.hg19.fasta \
--bam $cfdna_bam \
--variant-input $parents_vcf \
--only-use-input-alleles \
2>&1 \
>${output_vcf} \
| python /path/to/freebayes_hoobari_patch.py \
-b $cfdna_bam \
-parents_vcf $parents_vcf \
-m M12148W \
-p H12148W \
-r $region \
-db 12148

/path/to/bgzip -f $output_vcf
/path/to/tabix -f -p vcf ${output_vcf}.gz
