<h2> <img src="https://github.com/nshomron/hoobari/raw/master/misc/hoobari_logo.png" width=100/>, a cell-free DNA based fetal variant detector </h2>

## user manual and guide

--------

## Overview

*Hoobari* is a [Bayesian](http://en.wikipedia.org/wiki/Bayesian_inference) fetal variant caller designed to find SNPs (single-nucleotide polymorphisms) and indels (insertions and deletions) in a noninvasive manner.

*Hoobari* uses sequencing data from the mother, the father and most importantly from [cell-free DNA (cfDNA)](https://en.wikipedia.org/wiki/Cell-free_fetal_DNA) in the maternal plasma, which contains both fetal and maternal DNA fragments. Parental variant calling and pre-processing of the cfDNA are performed using [*FreeBayes*](https://github.com/ekg/freebayes).

## Obtaining

To download Hoobari:

    git clone --recursive git://github.com/nshomron/hoobari.git

Note the use of --recursive. This is required in order to download all nested git submodules for external repositories.

## Usage

Hoobari's pipeline consists of 3 steps:
1. Parental variant detection (Freebayes)
2. Pre-processing of cfDNA (Freebayes + patch)
3. Fetal variant calling (Hoobari)

**Parental variant detection:**
    
    freebayes \
    --fasta-reference h.sapiens.fasta \
    mother.sorted.mdup.bam \
    father.sorted.mdup.bam \
    | bgzip -c > parents.vcf.gz
    
    tabix -f -p vcf parents.vcf.gz

**Pre-processing of cfDNA:**
    
    #!/bin/bash

    freebayes \
    -d \
    --fasta-reference h.sapiens.fasta \
    --bam cfdna.sorted.mdup.bam \
    --variant-input parents.vcf.gz \
    --only-use-input-alleles \
    2>&1 \
    >cfdna.vcf \
    | python /path/to/hoobari/src/freebayes_patch.py \
    -b cfdna.sorted.mdup.bam \
    -parents_vcf parents.vcf.gz \
    -m MATERNAL_SAMPLE_NAME \
    -p PATERNAL_SAMPLE_NAME

    bgzip -f cfdna.vcf
    tabix -f -p vcf cfdna.vcf.gz

**Fetal variant calling:**
    
    hoobari \
    -m MATERNAL_SAMPLE_NAME \
    -p PATERNAL_SAMPLE_NAME \
    -f CFDNA_SAMPLE_NAME \
    -parents_vcf parents.vcf.gz \
    -cfdna_vcf cfdna.vcf.gz \
    > fetus.vcf
