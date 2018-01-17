#!/bin/bash
#$ -N freeb
#$ -S /bin/bash
#$ -cwd
#$ -o .$JOB_NAME.$JOB_ID.log

cat /groups/nshomron/guyshapira/tmp/hoobari/debug_trimmed.txt |python3 /groups/nshomron/guyshapira/projects/hoobari/src/freebayes_patch.py \
    -b /groups/nshomron/tomr/projects/cffdna/runs/fam01/short/S01.sorted.mdup.bam \
    -parents_vcf /groups/nshomron/guyshapira/tmp/hoobari/short.vcf.gz \
    -m M01 -p F01 -r chr1:0-1000000
