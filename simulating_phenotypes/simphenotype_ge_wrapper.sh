
#!/bin/bash

#script generates genetic effects and phenotypes with different patterns of stratification

if [ "$#" -ne 3 ]; then
    echo "usage: bash simphenotype_ge_wrapper.sh <path to pgen genotype file - prefix> <path to pop file> <path to output file>"
    exit 1
else

geno_file_prefix=${1} #full path of genotype file
pop_file=${2} # fill path of pop file
phenotype_output=${3} #path of phenotype output

#select 2000 variants to be causal

plink2 \
--bp-space 100000 \
--mac 1 \
--make-pgen \
--out ${geno_file_prefix}.thinned_100kb \
--pfile ${geno_file_prefix}

#calculate frequency of causal variants - this will be used to generate effects
plink2 --freq \
--out ${geno_file_prefix}.thinned_100kb.frq \
--pfile ${geno_file_prefix}.thinned_100kb

echo "simulating genetic effects"
#simulate effect sizes for these variants
Rscript simgeffects.R \
${geno_file_prefix}.thinned_100kb.frq.afreq \
${geno_file_prefix}.thinned_100kb.effects

#generate genetic values for each individual using their genotypes and these effects
plink2 --pfile ${geno_file_prefix}.thinned_100kb \
--out ${geno_file_prefix}.thinned_100kb.gvalue \
--score ${geno_file_prefix}.thinned_100kb.effects cols=dosagesum,scoresums

echo "simulating phenotypes"
#use these genetic values to generate phenotypes
Rscript simphenotype_ge.R \
${geno_file_prefix}.thinned_100kb.gvalue.sscore \
${pop_file} \
${phenotype_output}

fi
