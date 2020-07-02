
#!bin/bash

source /home/aazaidi/anaconda3/etc/profile.d/conda.sh
conda activate rpkgs

phenotype=${1}
iteration=${2}
mkdir -p train/ldsc.intercept/${phenotype}

echo "carryin out ldsc for cm pca"
Rscript ldsc_regression_noplots.R \
train/genos_complex_l1e7_ss500_m0.08_chr1_20.rmdup.train.ldsc.chr1_20.l2.ldscore \
train/gwas_results/fixed_effects/ge/gwas_complex_train.ge.${iteration}.cm.${phenotype}.glm.linear.gz \
train/genotypes/genos_complex_l1e7_ss500_m0.08_chr1_20.rmdup.train.snps.frq.afreq \
train/ldsc.intercept/${phenotype}/gwas_complex_train.ge.${iteration}.cm.${phenotype}

echo "carrying out ldsc for re pca"
Rscript ldsc_regression_noplots.R \
train/genos_complex_l1e7_ss500_m0.08_chr1_20.rmdup.train.ldsc.chr1_20.l2.ldscore \
train/gwas_results/fixed_effects/ge/gwas_complex_train.ge.${iteration}.re.${phenotype}.glm.linear.gz \
train/genotypes/genos_complex_l1e7_ss500_m0.08_chr1_20.rmdup.train.snps.frq.afreq \
train/ldsc.intercept/${phenotype}/gwas_complex_train.ge.${iteration}.re.${phenotype}

echo "carrying out ldsc for cmre pca"
Rscript ldsc_regression_noplots.R \
train/genos_complex_l1e7_ss500_m0.08_chr1_20.rmdup.train.ldsc.chr1_20.l2.ldscore \
train/gwas_results/fixed_effects/ge/gwas_complex_train.ge.${iteration}.cmre.${phenotype}.glm.linear.gz \
train/genotypes/genos_complex_l1e7_ss500_m0.08_chr1_20.rmdup.train.snps.frq.afreq \
train/ldsc.intercept/${phenotype}/gwas_complex_train.ge.${iteration}.cmre.${phenotype}
