<h2> <img src="https://github.com/nshomron/hoobari/raw/master/misc/hoobari_logo.png" width=100/>, a cell-free DNA based fetal variant detector </h2>

## user manual and guide

--------

## Overview

*Hoobari* is a [Bayesian](http://en.wikipedia.org/wiki/Bayesian_inference) genetic fetal variant detector designed to find small polymorphisms, specifically SNPs (single-nucleotide polymorphisms) and indels (insertions and deletions) smaller than the length of a short-read sequencing alignment.

*Hoobari* uses sequencing data from the mother, the father and most importantly from cell-free DNA (cfDNA) in the maternal plasma, which contains both fetal and maternal DNA fragments. Parental variant calling and pre-processing of the cfDNA are performed using [*FreeBayes*](https://github.com/ekg/freebayes).

## Obtaining

To download Hoobari:

    git clone --recursive git://github.com/nshomron/hoobari.git

Note the use of --recursive. This is required in order to download all nested git submodules for external repositories.

## Usage

Hoobari's pipeline consists of 3 steps:
1. Parental variant detection (Freebayes)
2. Pre-processing of cfDNA (Freebayes + patch)
3. Fetal variant calling (Hoobari)

To run Hoobari in the simplest way:

**Parental variant detection:**
    freebayes --fasta-reference h.sapiens.fasta mother.sorted.mdup.bam father.sorted.mdup.bam | bgzip -c > parents.vcf.gz
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
    -parents_vcf 02_parents.vcf.gz \
    -m M02 \
    -p F02 \
    -d

    bgzip -f cfdna.vcf
    tabix -f -p vcf cfdna.vcf.gz

**Fetal variant calling:**
    
    hoobari -m MATERNAL_SAMPLE_NAME -p PATERNAL_SAMPLE_NAME -f CFDNA_SAMPLE_NAME -parents_vcf parents.vcf.gz -cfdna_vcf cfdna.vcf.gz | bgzip -c > $out_vcf
    tabix -f -p vcf $out_vcf