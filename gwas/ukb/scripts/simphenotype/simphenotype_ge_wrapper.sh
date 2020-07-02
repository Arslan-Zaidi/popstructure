
#!/bin/bash

source /home/aazaidi/anaconda3/etc/profile.d/conda.sh
conda activate rpkgs

#script generates genetic effects and phenotypes with different patterns of stratification

if [ "$#" -lt 4 ]; then
    echo "usage: bash simphenotype_ge_wrapper.sh <path to pgen genotype file - prefix> <path to pop file> <path to output file> <replicate>"
    exit 1
else

geno_file_prefix=${1} #full path of genotype file
pop_file=${2} # fill path of pop file
phenotype_output=${3} #path of phenotype output
rep=${4} #replicate

#select 2000 variants to be causal

echo "simulating genetic effects"
#simulate effect sizes for these variants
Rscript scripts/simphenotype/simgeffects.R \
${geno_file_prefix}.snps.frq.afreq \
${geno_file_prefix}.${rep}.thinned_100kb.effects \
${rep}

#generate genetic values for each individual using their genotypes and these effects
plink2 --pfile ${geno_file_prefix} \
--out ${geno_file_prefix}.${rep}.thinned_100kb.gvalue \
--score ${geno_file_prefix}.${rep}.thinned_100kb.effects cols=dosagesum,scoresums

echo "simulating phenotypes"
#use these genetic values to generate phenotypes
Rscript scripts/simphenotype/simphenotype_ge.R \
${geno_file_prefix}.${rep}.thinned_100kb.gvalue.sscore \
${pop_file} \
${phenotype_output} \
${rep}

fi
