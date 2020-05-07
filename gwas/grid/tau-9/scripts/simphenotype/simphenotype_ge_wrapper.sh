
#!/bin/bash

source /home/aazaidi/anaconda3/etc/profile.d/conda.sh
conda activate rpkgs

#script generates genetic effects and phenotypes with different patterns of stratification

if [ "$#" -lt 3 ]; then
    echo "usage: bash simphenotype_ge_wrapper.sh <path to pgen genotype file - prefix> <path to pop file> <path to output file>"
    exit 1
else

geno_file_prefix=${1} #full path of genotype file
pop_file=${2} # fill path of pop file
phenotype_output=${3} #path of phenotype output
rep=${4}

#select 2000 variants to be causal

# plink2 \
# --bp-space 100000 \
# --seed ${rep} \
# --mac 1 \
# --make-pgen \
# --out ${geno_file_prefix}.${rep}.thinned_100kb \
# --pfile ${geno_file_prefix} \
# --read-freq ${geno_file_prefix}.frq.afreq

# #calculate frequency of causal variants - this will be used to generate effects
# plink2 --freq \
# --out ${geno_file_prefix}.${rep}.thinned_100kb.frq \
# --pfile ${geno_file_prefix}.${rep}.thinned_100kb

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

# #remove unnecessary files
# rm ${geno_file_prefix}.${rep}.thinned_100kb.pgen
# rm ${geno_file_prefix}.${rep}.thinned_100kb.psam
# rm ${geno_file_prefix}.${rep}.thinned_100kb.pvar

fi
