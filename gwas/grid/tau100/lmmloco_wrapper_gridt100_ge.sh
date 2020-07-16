#!/bin/bash

rep=${1}
chr=${2}
pheno_column=${3}

if [ $pheno_column -eq 1 ]; then
    phenotype="smooth"

elif [ $pheno_column -eq 3 ]; then
    phenotype="sharp"

else
  echo "pheno_ix can only be 1 (smooth) or 3 (sharp)"
  exit
fi

echo $pheno_column


mkdir -p train/gwas_results/mlma/ge

echo "running GWAS on heritable phenotypes"

#2. common GRM
bash gctaloco_mlma_gridt100_ge.sh \
cm.200k \
phenotypes/ge/pheno_gridt100_ge_s9k.train.${rep}.txt \
${pheno_column} \
${chr} \
train/gwas_results/mlma/ge/gwas_gridt100_train.ge.e${rep}.cm.${phenotype}.c${chr}

#3. rare GRM
bash gctaloco_mlma_gridt100_ge.sh \
re.all \
phenotypes/ge/pheno_gridt100_ge_s9k.train.${rep}.txt \
${pheno_column} \
${chr} \
train/gwas_results/mlma/ge/gwas_gridt100_train.ge.e${rep}.re.${phenotype}.c${chr}

echo "Zipping glm results"
#gzip -f train/gwas_results/mlma/ge/*.${rep}.*.mlma
