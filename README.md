<h1> h◎bari, a cell-free DNA based fetal variant detector </h1>

## user manual and guide

--------

## Overview

[*Hoobari*](https://github.com/nshomron/hoobari) is the first fetal variant calling program, designed to prenataly find SNPs (single-nucleotide polymorphisms) and indels (insertions and deletions) in a noninvasive manner. It requires sequencing data from the mother, the father and most importantly from [cell-free DNA (cfDNA)](https://en.wikipedia.org/wiki/Cell-free_fetal_DNA) in the maternal plasma, which contains both fetal and maternal DNA fragments.

*Hoobari* is based on a [Bayesian](http://en.wikipedia.org/wiki/Bayesian_inference) algorithm in which each cfDNA fragment has its own probability of being fetal. Its output is a standard Variant Calling Format ([VCF](https://samtools.github.io/hts-specs/VCFv4.2.pdf)) file, that can be further assessed by many bioinformatical tools.

One of *Hoobari*'s main goals is to create a general framework for the process of noninvasive fetal genotyping, that follows existing standards and is compatible with other bioinformatical tools. *Hoobari*'s workflow is therefore similar to that of other variant callers, and conveniently allows the introduction of future improvements. Hopefully, it will make this field more accessible for other researchers.

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
