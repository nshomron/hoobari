#!/bin/bash
#$ -N freeb
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -o .$JOB_NAME.$JOB_ID.log

python3 hoobari/src/freebayes_patch.py \
    -b /groups/nshomron/tomr/projects/cffdna/runs/fam01/short/S01.sorted.mdup.bam \
    -parents_vcf /groups/nshomron/guyshapira/tmp/hoobari/short.vcf \
    -m M01 -p F01 -r chr1:0-1000000
