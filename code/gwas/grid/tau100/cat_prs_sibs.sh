#!/bin/bash

pheno=${1}

echo "concatenating all"
for i in {1..20}; \
do tail -n+2 -q sibs_prs/${pheno}/gridt100_prs_p${pheno}.e${i}.all.sscore | \
awk -v OFS="\t" -v i=${i} '{print i,$0}'; done >sibs_prs/gridt100_prs_p${pheno}.all.sscore

echo "concatenating ascertained"
for i in {1..20}; \
do tail -n+2 -q sibs_prs/${pheno}/gridt100_prs_p${pheno}.e${i}.asc.nc.sscore | \
awk -v OFS="\t" -v i=${i} '{print i,$0}'; done >sibs_prs/gridt100_prs_p${pheno}.asc.nc.sscore
