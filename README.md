<h2> <img src="https://github.com/nshomron/hoobari/raw/master/misc/hoobari_logo.png" width=100/>, a cell-free DNA based fetal variant detector </h2>

## user manual and guide

--------

## Overview

*Hoobari* is a [Bayesian](http://en.wikipedia.org/wiki/Bayesian_inference) genetic fetal variant detector designed to find small polymorphisms, specifically SNPs (single-nucleotide polymorphisms) and indels (insertions and deletions) smaller than the length of a short-read sequencing alignment.

*Hoobari* uses sequencing data from the mother, the father and most importantly cell-free DNA (cfDNA) from the maternal plasma, which contains both fetal and maternal DNA fragments. Parental variant calling and pre-processing of the cfDNA are performed using [*FreeBayes*](https://github.com/ekg/freebayes).

## Obtaining

To download Hoobari:

    git clone --recursive git://github.com/nshomron/hoobari.git

Note the use of --recursive. This is required in order to download all nested git submodules for external repositories.

## Usage

Hoobari's pipeline consists of 3 simple steps:
1. Parental variant detection (Freebayes)
2. Pre-processing of cfDNA (Freebayes + patch)
3. Fetal variant calling (Hoobari)

