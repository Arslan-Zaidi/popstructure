#!/bin/bash

chrom=${1}

#output imputation accuracy scores
zcat train/imputed_genotypes/imputed_chr${chrom}.beagle.vcf.gz | \
grep -v '^#' | cut -f3,8 > train/imputed_genotypes/imputed_chr${chrom}.beagle.dr2
