#!/bin/bash

## batch-queuing system arguments
#$ -N jobname
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -o .$JOB_NAME.$JOB_ID.log

# terminate on error
set -e

region=$1
vcf_dir='M12148_plasma'
output_vcf=${vcf_dir}/M12148_plasma.${region}.vcf
mkdir -p ${vcf_dir}
mkdir -p tmp_hb

. /share/apps/modules/init/bash

module load python/anaconda3-4.0.0

freebayes \
-d \
--region $region \
--fasta-reference ~/nshomron_tomr/tools/GATK/bundle/ucsc/ucsc.hg19.fasta \
--bam M12148_plasma.sorted.mdup.bam \
--variant-input 12148_parents_wgs_all.vcf.gz \
--only-use-input-alleles \
2>&1 \
>${output_vcf} \
| python freebayes_dd_hoobari_patch.py \
-b M12148_plasma.sorted.mdup.bam \
-parents_vcf 12148_parents_wgs_all.vcf.gz \
-m M12148W \
-p H12148W \
-r $region \
-db 12148

bgzip -f $output_vcf
tabix -f -p vcf ${output_vcf}.gz
