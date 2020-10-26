
#!/bin/bash

source /home/aazaidi/anaconda3/etc/profile.d/conda.sh
conda activate rpkgs

#script generates non-genetic effects and phenotypes with different patterns of stratification

if [ "$#" -lt 3 ]; then
    echo "usage: bash simphenotype_noge_wrapper.sh <path to pop file> <path to output file> <replicate>"
    exit 1
else

#geno_file_prefix=${1} #full path of genotype file
pop_file=${1} # fill path of pop file
phenotype_output=${2} #path of phenotype output
rep=${3} #replicate

echo "simulating phenotypes"
Rscript scripts/simphenotype/simphenotype_noge.R \
${pop_file} \
${phenotype_output} \
${rep}

fi
