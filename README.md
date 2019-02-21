<h1> h◎bari </h1>
<h2> a Bayesian-based noninvasive fetal variant detector </h2>

## user manual and guide

--------

## Overview

[*Hoobari*](https://github.com/nshomron/hoobari) *(Hebrew: hoo - him, bari - healthy, oobari - fetal)* is the first fetal variant calling program, designed to prenataly find SNPs (single-nucleotide polymorphisms) and indels (insertions and deletions) in a noninvasive manner. It requires sequencing data of the mother, the father and the [cell-free DNA (cfDNA)](https://en.wikipedia.org/wiki/Cell-free_fetal_DNA), that is found in the maternal plasma and contains both fetal and maternal DNA fragments.

*Hoobari* is based on a [Bayesian](http://en.wikipedia.org/wiki/Bayesian_inference) algorithm in which each cfDNA fragment has its own probability of being fetal. Its output is a standard Variant Calling Format ([VCF](https://samtools.github.io/hts-specs/VCFv4.2.pdf)) file, that can be further analyzed and annonated by different tools.

One of *Hoobari*'s main goals is to create a general framework for the process of noninvasive fetal genotyping, that follows existing standards and is compatible with other bioinformatical tools. *Hoobari*'s workflow is therefore similar to that of other variant callers, and conveniently allows the introduction of future improvements. Hopefully, *Hoobari* will help make this field, which is still in its infancy, somewhat more accessible for other researchers.

## Citing Hoobari

If you are using the *Hoobari* in your research, please cite our paper as follows:

Bayesian-based noninvasive prenatal diagnosis of single-gene disorders. Tom Rabinowitz, Avital Polsky, David Golan, Artem Danilevsky, Guy Shapira, Chen Raff, Lina Basel-Salmon, Reut Tomashov Matar, and Noam Shomron. *Genome Research*. 2019. [doi:10.1101/gr.235796.118](https://genome.cshlp.org/content/early/2019/02/13/gr.235796.118.abstract)

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
    > parents.vcf

**Pre-processing of cfDNA:**
    
    freebayes \
    -d \
    --fasta-reference h.sapiens.fasta \
    --bam cfdna.sorted.mdup.bam \
    --variant-input parents.vcf \
    --only-use-input-alleles \
    |& python /path/to/hoobari/src/freebayes_patch.py \
    -b cfdna.sorted.mdup.bam \
    -parents_vcf parents.vcf \
    -m mother.sorted.mdup.bam \
    -p father.sorted.mdup.bam

**Fetal variant calling:**
    
    hoobari \
    -parents_vcf parents.vcf \
    -cfdna_vcf cfdna.vcf
